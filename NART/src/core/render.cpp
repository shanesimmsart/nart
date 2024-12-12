#include <atomic>
#include <thread>

#include "render.h"

#define FilterTableResolution 64

#define DEBUG_BUCKET 0
#define DEBUG_BUCKET_X 36
#define DEBUG_BUCKET_Y 22

#define BSDF_SAMPLING 1
#define LIGHT_SAMPLING 1
#define MAX_BOUNCES 10

RenderSession::RenderSession(const Scene& scene, uint32_t imageWidth, uint32_t imageHeight, uint32_t bucketSize, uint32_t spp, float filterWidth, float rougheningFactor) :
    scene(scene), imageWidth(imageWidth), imageHeight(imageHeight), bucketSize(bucketSize), spp(spp), filterWidth(filterWidth), rougheningFactor(rougheningFactor)
{
    filterBounds = static_cast<uint32_t>(glm::ceil(filterWidth));
    tileSize = bucketSize + (filterBounds * 2);
    totalWidth = imageWidth + (filterBounds * 2);
    totalHeight = imageHeight + (filterBounds * 2);
}

void RenderSession::AddSample(const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels)
{
    // Calculate discrete x and y bounds of filter
    uint32_t x0 = static_cast<uint32_t>(glm::floor(sampleCoords.x - filterWidth));
    uint32_t x1 = static_cast<uint32_t>(glm::ceil(sampleCoords.x + filterWidth));
    uint32_t y0 = static_cast<uint32_t>(glm::floor(sampleCoords.y - filterWidth));
    uint32_t y1 = static_cast<uint32_t>(glm::ceil(sampleCoords.y + filterWidth));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Calculate distance from sampleCoords
            float distX = (static_cast<float>(x) + 0.5f) - sampleCoords.x;
            float distY = (static_cast<float>(y) + 0.5f) - sampleCoords.y;
            float dist = glm::sqrt(distX * distX + distY * distY);

            // Get index into precomputed filter table based on distance from sample, as opposed to:
            // float filterWeight = Gaussian(filterWidth, dist);
            uint8_t filterIndex = glm::min(static_cast<uint8_t>((dist / filterWidth) * FilterTableResolution), static_cast<uint8_t>(FilterTableResolution - 1));
            float filterWeight = filterTable[filterIndex];

            // Transform image coords into discrete tile coords
            uint32_t tileX = static_cast<uint32_t>(glm::floor(distX + glm::mod(sampleCoords.x - static_cast<float>(filterBounds), static_cast<float>(bucketSize)) + static_cast<float>(filterBounds)));
            uint32_t tileY = static_cast<uint32_t>(glm::floor(distY + glm::mod(sampleCoords.y - static_cast<float>(filterBounds), static_cast<float>(bucketSize)) + static_cast<float>(filterBounds)));

            uint32_t tileIndex = (tileY * tileSize) + tileX;

            pixels[tileIndex].contribution += L * filterWeight;
            pixels[tileIndex].filterWeightSum += filterWeight;
        }
    }
}

glm::vec3 RenderSession::EstimateDirect(const glm::vec3 wo, BSDF& bsdf, const Intersection& isect, float* alphaTweak, const Ray& ray, RNG& rng, uint8_t* flags)
{
    glm::vec3 L = glm::vec3(0.f);

    glm::vec3 wi;
    glm::vec3 Li;
    float scatteringPdf = 0.f;
    float lightingPdf = 0.f;

    float numLights = scene.GetNumLights();
    uint8_t lightIndex = static_cast<uint8_t>(glm::min(rng.UniformFloat(), 1.f - glm::epsilon<float>()) * numLights);

    const Light& light = scene.GetLight(lightIndex);
    glm::vec3 f(0.f);

    // Compute direct light
#if BSDF_SAMPLING
    scatteringPdf = 0.f;
    glm::vec2 scatterSample(rng.UniformFloat(), rng.UniformFloat());
    float bsdfSample = rng.UniformFloat();

    f = bsdf.Sample_f(wo, &wi, bsdfSample, scatterSample, &scatteringPdf, flags, 1);

    if (scatteringPdf > 0.f)
    {
        float flip = wi.z > 0.f ? 1.f : -1.f;
        Ray shadowRay(isect.p + (isect.gn * shadowBias * flip), bsdf.ToWorld(wi));
        Intersection lightIsect;
        Li = light.Li(&lightIsect, isect.p, bsdf.ToWorld(wi), &lightingPdf);
        if (!scene.GetBVH().Intersect(shadowRay, &lightIsect))
        {
            float weight = 1.f;

            if (!(*flags & SPECULAR))
            {
#if LIGHT_SAMPLING
                weight = (scatteringPdf * scatteringPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                if (lightingPdf > 0.f)
                {
                    L += (f * Li * glm::abs(wi.z) * weight) / scatteringPdf;
                }
            }

            else
            {
                L += (f * Li * glm::abs(wi.z) * weight) / scatteringPdf;
            }
        }
    }
#endif
#if LIGHT_SAMPLING
    glm::vec3 wiWorld;
    lightingPdf = 0.f;
    glm::vec2 lightSample(rng.UniformFloat(), rng.UniformFloat());

    // lightIsect.tMax used for checking shadows
    Intersection lightIsect;
    Li = light.Sample_Li(&lightIsect, isect.p, &wiWorld, lightSample, &lightingPdf);
    wi = bsdf.ToLocal(wiWorld);

    float flip = wi.z > 0.f ? 1.f : -1.f;
    Ray shadowRay(isect.p + (isect.gn * shadowBias * flip), wiWorld);
    if (!scene.GetBVH().Intersect(shadowRay, &lightIsect) && lightingPdf > 0.f) // || flags && TRANSMISSIVE)
    {
        float weight = 1.f;
        scatteringPdf = bsdf.Pdf(wo, wi, 1);
        if (scatteringPdf > 0.f)
        {
            glm::vec3 f = bsdf.f(wo, wi, 1);
#if BSDF_SAMPLING
            weight = (lightingPdf * lightingPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
            L += (f * Li * glm::abs(wi.z) * weight) / lightingPdf;
        }
    }
#endif

    return L * numLights;
}

std::vector<Pixel> RenderSession::RenderTile(const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1)
{
    std::vector<Pixel> pixels((tileSize) * (tileSize));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Create stratified image samples in a Latin square pattern
            RNG rng;
            rng.Seed(y * totalWidth + x);

            std::vector<glm::vec2> imageSamples(spp);
            LatinSquare(rng, spp, imageSamples);

            for (uint32_t i = 0; i < spp; ++i)
            {
                glm::vec2 imageSample = imageSamples[i];

                const Camera& camera = scene.GetCamera();
                Ray ray = camera.CastRay(imageSample, imageWidth, imageHeight, x, y);

                // Radiance
                glm::vec3 L(0.f, 0.f, 0.f);
                float alpha = 0.f;
                Intersection isect;
                // Throughput
                glm::vec3 beta(1.f);
                uint8_t flags;

                float gamma = rougheningFactor * rougheningFactor;
                float alphaTweak = 1.f;

                for (uint32_t bounce = 0; bounce < MAX_BOUNCES; ++bounce)
                {
                    // Light intersections
                    float lightTMax = isect.tMax;
                    bool lightHit = false;
                    glm::vec3 Le(0.f);
                    for (uint8_t i = 0; i < scene.GetNumLights(); ++i)
                    {
                        Intersection lightIsect;
                        const Light& light = scene.GetLight(i);
                        glm::vec3 Li = light.Li(&lightIsect, ray.o, ray.d);
                        if (lightIsect.tMax < lightTMax)
                        {
                            Le = Li;
                            lightTMax = lightIsect.tMax;
                            isect.tMax = lightIsect.tMax;
                            lightHit = true;
                            alpha = 1.f;
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

                        BSDF bsdf = isect.material->CreateBSDF(isect.sn, alphaTweak);
                        glm::vec3 wo = bsdf.ToLocal(-ray.d);

                        L += EstimateDirect(wo, bsdf, isect, &alphaTweak, ray, rng, &flags) * beta;

                        // Spawn new ray
                        glm::vec2 scatteringSample(rng.UniformFloat(), rng.UniformFloat());
                        float bsdfSample = rng.UniformFloat();

                        float scatteringPdf = 0.f;
                        // Sample for new alpha
                        float alpha_i;
                        glm::vec3 wi;
                        glm::vec3 f = bsdf.Sample_f(wo, &wi, bsdfSample, scatteringSample, &scatteringPdf, &flags, 0, &alpha_i);
                        // Sample for new ray
                        if (scatteringPdf <= 0.f) break;
                        alphaTweak = (1.f - (gamma * alpha_i)) * alphaTweak;
                        beta *= (f / scatteringPdf) * glm::abs(wi.z);
                        // Transform to world
                        float flip = wi.z > 0.f ? 1.f : -1.f;
                        ray = Ray(isect.p + (isect.gn * shadowBias * flip), bsdf.ToWorld(wi));

                        // Russian roulette
                        float q = glm::max((beta.x + beta.y + beta.z) * 0.33333f, 0.f);

                        if (bounce > 3)
                        {
                            if (q >= rng.UniformFloat())
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
                glm::vec2 sampleCoords = glm::vec2(static_cast<float>(x + filterBounds) + imageSample.x, static_cast<float>(y + filterBounds) + imageSample.y);
                // Add sample to pixels within filter width
                AddSample(filterTable, sampleCoords, glm::vec4 (L, alpha), pixels);
            }
        }
    }

    return pixels;
}

std::vector<Pixel> RenderSession::Render()
{
    std::vector<Pixel> pixels(totalHeight * totalWidth);

    uint32_t nBucketsX = uint32_t(glm::ceil(static_cast<float>(imageWidth) / static_cast<float>(bucketSize)));
    uint32_t nBucketsY = uint32_t(glm::ceil(static_cast<float>(imageHeight) / static_cast<float>(bucketSize)));
    uint32_t nBuckets = nBucketsX * nBucketsY;

    // Pre-compute filter values (only need to do this in 1D as filter is currently only isotropic)
    float filterTable[FilterTableResolution];
    for (uint8_t i = 0; i < FilterTableResolution; ++i)
    {
        filterTable[i] = Gaussian(FilterTableResolution - 1, i);
    }

    // Create tiles and render each one in parallel
    std::vector<Pixel> p(tileSize * tileSize);
    std::vector<std::vector<Pixel>> tiles(nBuckets, p);

    std::atomic<uint32_t> nBucketsComplete = 0;
#if !DEBUG_BUCKET
    std::thread progressLogger([&nBucketsComplete, nBuckets]
        {
            while (nBucketsComplete < nBuckets)
            {
                std::cout << "\r" << static_cast<int>(glm::floor((static_cast<float>(nBucketsComplete) / static_cast<float>(nBuckets)) * 100.f)) << "%" << std::flush;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            std::cout << "\r100%" << std::endl; // Ensure final output is 100%
    });
#endif

    tbb::task_group tg;

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
            tg.run([this, &filterTable, &tiles, index, x, y, &nBucketsComplete, nBuckets]
                {
                    uint32_t x0 = bucketSize * x;
                    uint32_t y0 = bucketSize * y;
                    uint32_t x1 = glm::min(bucketSize * (x + 1), totalWidth);
                    uint32_t y1 = glm::min(bucketSize * (y + 1), totalHeight);
                    std::vector<Pixel> tile = RenderTile(filterTable, x0, x1, y0, y1);
                    tiles[index] = tile;
                    nBucketsComplete++;
                });
        }
    }

    tg.wait();
#if !DEBUG_BUCKET
    progressLogger.join();
#endif

    // Combine tiles into image
    for (uint32_t j = 0; j < nBucketsY; ++j)
    {
        for (uint32_t i = 0; i < nBucketsX; ++i)
        {
            for (uint32_t y = 0; y < tileSize; ++y)
            {
                for (uint32_t x = 0; x < tileSize; ++x)
                {
                    std::vector<Pixel> v = tiles[j * nBucketsX + i];
                    uint32_t pX = x + (i * bucketSize);
                    uint32_t pY = y + (j * bucketSize);

                    // Ignore tile pixels beyond image bounds
                    if (pX < (imageWidth + filterBounds) && pY < (imageHeight + filterBounds))
                    {
                        uint32_t pIndex = (pY * totalWidth) + pX;
                        pixels[pIndex].contribution += v[(y * tileSize) + x].contribution;
                        pixels[pIndex].filterWeightSum += v[(y * tileSize) + x].filterWeightSum;
                    }
                }
            }
        }
    }

    return pixels;
}

void RenderSession::WriteImageToEXR(const std::vector<Pixel>& pixels, const char* filePath)
{
    // Write image to EXR
    Imf::Array2D<Imf::Rgba> imagePixels(imageHeight, imageWidth);

    // Ignore pixels beyond image bounds
    for (uint32_t y = filterBounds; y < (imageHeight + filterBounds); ++y)
    {
        for (uint32_t x = filterBounds; x < (imageWidth + filterBounds); ++x)
        {
            glm::vec4 result(0.f);

            result = pixels[y * totalWidth + x].contribution / pixels[y * totalWidth + x].filterWeightSum;

            imagePixels[y - filterBounds][x - filterBounds].r = result.r;
            imagePixels[y - filterBounds][x - filterBounds].g = result.g;
            imagePixels[y - filterBounds][x - filterBounds].b = result.b;
            imagePixels[y - filterBounds][x - filterBounds].a = result.a;
        }
    }

    Imf::RgbaOutputFile file(filePath, imageWidth, imageHeight, Imf::WRITE_RGBA);
    file.setFrameBuffer(&imagePixels[0][0], 1, imageWidth);
    file.writePixels(imageHeight);
}

std::vector<std::unique_ptr<RenderSession>> LoadSessions(std::string scenePath, const Scene& scene)
{
    nlohmann::json json;
    std::ifstream ifs;
    ifs.open(scenePath);

    if (!ifs)
    {
        std::cerr << "Error: Scene file could not be opened.\n";
    }

    try
    {
        ifs >> json;
    }

    catch (nlohmann::json::exception& e)
    {
        std::cerr << "Error parsing scene: " << e.what() << "\nAborting.\n";
        std::abort();
    }

    std::vector<std::unique_ptr<RenderSession>> sessions;

    if (!json["renderSessions"].is_null()) {
        for (auto& elem : json["renderSessions"])
        {
            // Default values if not found in JSON
            uint32_t imageWidth = 64;
            uint32_t imageHeight = 64;
            uint32_t bucketSize = 32;
            uint32_t spp = 1;
            float filterWidth = 1.f;
            float rougheningFactor = 0.f;

            if (!elem["imageWidth"].is_null()) imageWidth = elem["imageWidth"].get<uint32_t>();
            if (!elem["imageHeight"].is_null()) imageHeight = elem["imageHeight"].get<uint32_t>();
            if (!elem["bucketSize"].is_null()) bucketSize = elem["bucketSize"].get<uint32_t>();
            if (!elem["spp"].is_null()) spp = elem["spp"].get<uint32_t>();
            if (!elem["filterWidth"].is_null()) filterWidth = elem["filterWidth"].get<float>();
            if (!elem["rougheningFactor"].is_null()) rougheningFactor = glm::min(glm::max(elem["rougheningFactor"].get<float>(), 0.f), 1.f);

            sessions.push_back(std::make_unique<RenderSession>(scene, imageWidth, imageHeight, bucketSize, spp, filterWidth, rougheningFactor));
        }
    }

    return sessions;
}


