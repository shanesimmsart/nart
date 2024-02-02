#include "render.h"

#define FilterTableResolution 64

#define DEBUG_BUCKET 0
#define DEBUG_BUCKET_X 53
#define DEBUG_BUCKET_Y 21

#define BSDF_SAMPLING 1
#define LIGHT_SAMPLING 1

float Gaussian(float width, float x)
{
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * glm::pi<float>() * sigma * sigma)) * glm::exp(-(x * x) / (2.f * sigma * sigma));
}

void AddSample(const RenderInfo& info, const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels)
{
    // Calculate discrete x and y bounds of filter
    uint32_t x0 = static_cast<uint32_t>(glm::floor(sampleCoords.x - info.filterWidth));
    uint32_t x1 = static_cast<uint32_t>(glm::ceil(sampleCoords.x + info.filterWidth));
    uint32_t y0 = static_cast<uint32_t>(glm::floor(sampleCoords.y - info.filterWidth));
    uint32_t y1 = static_cast<uint32_t>(glm::ceil(sampleCoords.y + info.filterWidth));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Calculate distance from sampleCoords
            float distX = (static_cast<float>(x) + 0.5f) - sampleCoords.x;
            float distY = (static_cast<float>(y) + 0.5f) - sampleCoords.y;
            float dist = glm::sqrt(distX * distX + distY * distY);

            // Get index into precomputed filter table based on distance from sample, as opposed to:
            // float filterWeight = Gaussian(scene.info.filterWidth, dist);
            uint8_t filterIndex = glm::min(static_cast<uint8_t>((dist / info.filterWidth) * FilterTableResolution), static_cast<uint8_t>(FilterTableResolution - 1));
            float filterWeight = filterTable[filterIndex];

            // Transform image coords into discrete tile coords
            uint32_t tileX = static_cast<uint32_t>(glm::floor(distX + glm::mod(sampleCoords.x - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));
            uint32_t tileY = static_cast<uint32_t>(glm::floor(distY + glm::mod(sampleCoords.y - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));

            uint32_t tileIndex = (tileY * info.tileSize) + tileX;

            pixels[tileIndex].contribution += L * filterWeight;
            pixels[tileIndex].filterWeightSum += filterWeight;
        }
    }
}

std::vector<Pixel> RenderTile(const Scene& scene, const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1)
{
    std::vector<Pixel> pixels((scene.info.tileSize) * (scene.info.tileSize));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Don't render beyond total image extents when the current tile goes beyond them
            if (x < scene.info.totalWidth && y < scene.info.totalHeight)
            {
                // Create stratified image samples in a Latin square pattern
                std::default_random_engine rng;
                rng.seed(y * scene.info.totalWidth + x);
                std::vector<glm::vec2> imageSamples(scene.info.spp);
                LatinSquare(rng, scene.info.spp, imageSamples);

                for (uint32_t i = 0; i < scene.info.spp; ++i)
                {
                    glm::vec2 imageSample = imageSamples[i];

                    // Create ray in camera-space
                    float aspectRatio = static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.imageHeight);
                    float px = (((static_cast<float>(x) + imageSample.x) / static_cast<float>(scene.info.imageWidth)) * 2.f - 1.f) * glm::tan(glm::radians(scene.camera.fov)) * aspectRatio;
                    float py = (((static_cast<float>(y) + imageSample.y) / static_cast<float>(scene.info.imageHeight)) * -2.f + 1.f) * glm::tan(glm::radians(scene.camera.fov));

                    glm::vec4 o = glm::vec4(0.f, 0.f, 0.f, 1.f);
                    glm::vec4 d(px, py, -1.f, 0.f);
                    glm::normalize(d);

                    // Transform ray into world-space using camera transform
                    glm::mat4 camToWorld = scene.camera.cameraToWorld;

                    o = o * camToWorld;
                    d = d * camToWorld;

                    Ray ray(o, d);

                    // Radiance
                    glm::vec3 L(0.f, 0.f, 0.f);
                    float alpha = 0.f;
                    Intersection isect;
                    // Throughput
                    glm::vec3 beta(1.f);
                    BSDFFlags flags;

                    float roughnessOffset = 0.f;

                    for (uint32_t bounce = 0; bounce < 10; ++bounce)
                    {
                        // TODO: Move this inside of scene.Intersect()
                        float lightTMax = isect.tMax;
                        bool lightHit = false;
                        glm::vec3 Le(0.f);
                        for (const auto& light : scene.lights)
                        {
                            Intersection lightIsect;
                            glm::vec3 Li = light->Li(&lightIsect, ray.o, ray.d);
                            if (lightIsect.tMax < lightTMax)
                            {
                                Le = Li;
                                lightTMax = lightIsect.tMax;
                                isect.tMax = lightIsect.tMax;
                                lightHit = true;
                            }
                        }

                        if (bounce == 0)
                        {
                            if (lightHit == true)
                            {
                                L = Le;
                                break;
                            }
                        }

                        if (scene.Intersect(ray, &isect))
                        {
                            if (bounce == 0) alpha = 1.f;

                            std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());

                            BSDF bsdf = isect.material->CreateBSDF(isect.sn, roughnessOffset);
                            glm::vec3 wo = bsdf.ToLocal(-ray.d);

                            float shadowBias = 0.0001f;

                            glm::vec3 wi;
                            glm::vec3 Li;
                            float scatteringPdf = 0.f;
                            float lightingPdf = 0.f;

                            float numLights = static_cast<float>(scene.lights.size());
                            uint8_t lightIndex = static_cast<uint8_t>(glm::min(distribution(rng), 1.f - glm::epsilon<float>()) * numLights);
                            std::shared_ptr<Light> light = scene.lights[lightIndex];
                            glm::vec3 f(0.f);

                            // Compute direct light
#if BSDF_SAMPLING
                            scatteringPdf = 0.f;
                            glm::vec2 scatterSample(distribution(rng), distribution(rng));
                            f = bsdf.Sample_f(wo, &wi, scatterSample, &scatteringPdf, &flags);

                            if (scatteringPdf > 0.f)
                            {
                                Ray shadowRay(isect.p + (isect.gn * shadowBias), bsdf.ToWorld(wi));
                                Intersection lightIsect;
                                Li = light->Li(&lightIsect, isect.p, bsdf.ToWorld(wi));
                                if (!scene.bvh->Intersect(shadowRay, &lightIsect))
                                {
                                    float weight = 1.f;

                                    if (!(flags & SPECULAR))
                                    {
                                        // TODO: I should probably do this in one function
                                        lightingPdf = light->Pdf(&lightIsect, isect.p, bsdf.ToWorld(wi));
#if LIGHT_SAMPLING
                                        weight = (scatteringPdf * scatteringPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                                        if (lightingPdf > 0.f)
                                        {
                                            L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / scatteringPdf;
                                            L *= numLights;
                                        }
                                    }

                                    else
                                    {
                                        L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / scatteringPdf;
                                        L *= numLights;
                                    }
                                }
                            }
#endif
#if LIGHT_SAMPLING
                            glm::vec3 wiWorld;
                            lightingPdf = 0.f;
                            glm::vec2 lightSample(distribution(rng), distribution(rng));
                            // lightIsect.tMax used for checking shadows
                            Intersection lightIsect;
                            Li = light->Sample_Li(&lightIsect, isect.p, &wiWorld, lightSample, &lightingPdf);

                            Ray shadowRay(isect.p + (isect.gn * shadowBias), wiWorld);
                            if (!scene.bvh->Intersect(shadowRay, &lightIsect) && lightingPdf > 0.f)
                            {
                                float weight = 1.f;
                                wi = bsdf.ToLocal(wiWorld);
                                scatteringPdf = bsdf.Pdf(wo, wi);
                                if (scatteringPdf > 0.f)
                                {
                                    glm::vec3 f = bsdf.f(wo, wi);
#if BSDF_SAMPLING
                                    weight = (lightingPdf * lightingPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                                    L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / lightingPdf;
                                    // L *= numLights;
                                }
                            }
#endif

                            // Spawn new ray
                            glm::vec2 scatteringSample(distribution(rng), distribution(rng));
                            scatteringPdf = 0.f;
                            f = bsdf.Sample_f(wo, &wi, scatteringSample, &scatteringPdf, &flags);
                            if (scatteringPdf <= 0.f) break;
                            if (flags & DIFFUSE) roughnessOffset += 0.5f;
                            beta *= (f / scatteringPdf) * glm::abs(wi.z);
                            // Transform to world
                            ray = Ray(isect.p + (isect.gn * shadowBias), bsdf.ToWorld(wi));

                            // Russian roulette
                            float q = glm::max((beta.x + beta.y + beta.z) * 0.33333f, 0.f);

                            // std::cout << q << "\n";

                            if (bounce > 3)
                            {
                                if (q >= distribution(rng))
                                {
                                    beta /= q;
                                }

                                else break;
                            }

                            isect = Intersection();
                        }

                        else break;
                    }

                    // Transform image sample to "total" image coords (image including filter bounds)
                    glm::vec2 sampleCoords = glm::vec2(static_cast<float>(x + scene.info.filterBounds) + imageSample.x, static_cast<float>(y + scene.info.filterBounds) + imageSample.y);
                    // Add sample to pixels within filter width
                    // Sample coords are respective to total image size
                    glm::vec4 Lalpha(L, alpha);
                    AddSample(scene.info, filterTable, sampleCoords, Lalpha, pixels);
                }
            }
        }
    }

    return pixels;
}

std::vector<Pixel> Render(const Scene& scene)
{
    std::vector<Pixel> pixels(scene.info.totalHeight * scene.info.totalWidth);

    uint32_t nBucketsX = uint32_t(glm::ceil(static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBucketsY = uint32_t(glm::ceil(static_cast<float>(scene.info.imageHeight) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBuckets = nBucketsX * nBucketsY;

    // Pre-compute filter values (only need to do this in 1D as filter is currently only isotropic)
    float filterTable[FilterTableResolution];
    for (uint8_t i = 0; i < FilterTableResolution; ++i)
    {
        filterTable[i] = Gaussian(FilterTableResolution - 1, i);
    }

    // Create tiles and render each one in parallel
    std::vector<Pixel> p(scene.info.tileSize * scene.info.tileSize);
    std::vector<std::vector<Pixel>> tiles(nBuckets, p);

    tbb::task_group tg;

    uint32_t nBucketsComplete = 0;

#if DEBUG_BUCKET
    for (uint32_t y = DEBUG_BUCKET_Y; y < DEBUG_BUCKET_Y + 1; ++y)
    {
        for (uint32_t x = DEBUG_BUCKET_X; x < DEBUG_BUCKET_X + 1; ++x)
        {
#else
    for (uint32_t y = 0; y < nBucketsY; ++y)
    {
        for (uint32_t x = 0; x < nBucketsX; ++x)
        {
#endif
            uint32_t index = y * nBucketsX + x;
            tg.run([&scene, &filterTable, &tiles, index, x, y, &nBucketsComplete, nBuckets]
                {
                    std::vector<Pixel> tile = RenderTile(scene, filterTable, scene.info.bucketSize * x, scene.info.bucketSize * (x + 1), scene.info.bucketSize * y, scene.info.bucketSize * (y + 1));
                    tiles[index] = tile;
                    nBucketsComplete++;
                    // TODO: Multi-threaded logging?
                    std::cout << "\r" << (int)glm::floor(((float)nBucketsComplete / (float)nBuckets) * 100.f) << "%   " << std::flush;
                });
        }
    }

    tg.wait();

    // Combine tiles into image
    std::cout << "\nCombining tiles into image...\n";

    for (uint32_t j = 0; j < nBucketsY; ++j)
    {
        for (uint32_t i = 0; i < nBucketsX; ++i)
        {
            for (uint32_t y = 0; y < scene.info.tileSize; ++y)
            {
                for (uint32_t x = 0; x < scene.info.tileSize; ++x)
                {
                    std::vector<Pixel> v = tiles[j * nBucketsX + i];
                    uint32_t pX = x + (i * scene.info.bucketSize);
                    uint32_t pY = y + (j * scene.info.bucketSize);
                    // Ignore tile pixels beyond image bounds
                    if (pX < (scene.info.imageWidth + scene.info.filterBounds) && pY < (scene.info.imageHeight + scene.info.filterBounds))
                    {
                        uint32_t pIndex = (pY * scene.info.totalWidth) + pX;
                        pixels[pIndex].contribution += v[(y * scene.info.tileSize) + x].contribution;
                        pixels[pIndex].filterWeightSum += v[(y * scene.info.tileSize) + x].filterWeightSum;
                    }
                }
            }
        }
    }

    return pixels;
}

void WriteImageToEXR(const RenderInfo& info, const std::vector<Pixel>& pixels, const char* filePath)
{
    // Write image to EXR
    Imf::Array2D<Imf::Rgba> imagePixels(info.imageHeight, info.imageWidth);

    // Ignore pixels beyond image bounds
    for (uint32_t y = info.filterBounds; y < (info.imageHeight + info.filterBounds); ++y)
    {
        for (uint32_t x = info.filterBounds; x < (info.imageWidth + info.filterBounds); ++x)
        {
            glm::vec4 result(0.f);

            result = pixels[y * info.totalWidth + x].contribution / pixels[y * info.totalWidth + x].filterWeightSum;

            imagePixels[y - info.filterBounds][x - info.filterBounds].r = result.r;
            imagePixels[y - info.filterBounds][x - info.filterBounds].g = result.g;
            imagePixels[y - info.filterBounds][x - info.filterBounds].b = result.b;
            imagePixels[y - info.filterBounds][x - info.filterBounds].a = result.a;
        }
    }

    Imf::RgbaOutputFile file(filePath, info.imageWidth, info.imageHeight, Imf::WRITE_RGBA);
    file.setFrameBuffer(&imagePixels[0][0], 1, info.imageWidth);
    file.writePixels(info.imageHeight);
}


