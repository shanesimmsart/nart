#include "../../include/nart/integrators/pathintegrator.h"

#define BSDF_SAMPLING 1
#define LIGHT_SAMPLING 1
#define MAX_BOUNCES 10

bool PathIntegrator::IsectIsValid(
    const Intersection& isect,
    const std::vector<IntersectionInfo, ArenaAllocator<IntersectionInfo>>&
        isectList,
    float& eta_outer) const {
    // Update outer eta for nested dielectrics
    eta_outer = 1.f;

    if (!isectList.empty())
        // If the medium we previously entered isn't the same as
        // the one we just intersected, use that mediums IOR
        if (isectList.end()[-1].meshID != isect.meshID)
            eta_outer = isectList.end()[-1].eta;
        // If it is the same we are exiting that medium,
        // so set the IOR to the penultimate medium in
        // the list, or keep IOR at 1 if there's nothing
        // else in the list
        else if (isectList.size() >= 2)
            eta_outer = isectList.end()[-2].eta;

    // Check if intersection is lower priority and
    // should therefore be ignored for nested dielectrics
    for (auto elem : isectList) {
        if (isect.priority < elem.priority) {
            return false;
        }
    }

    return true;
}

glm::vec3 PathIntegrator::EstimateDirect(const Scene& scene, const glm::vec3 wo,
                                         BSDF& bsdf, const Intersection& isect,
                                         const Ray& ray, RNG& rng,
                                         uint8_t& flags,
                                         float eta_outer) const {
    glm::vec3 L = glm::vec3(0.f);

    glm::vec3 wi;
    glm::vec3 Li;
    float scatteringPdf = 0.f;
    float lightingPdf = 0.f;

    float numLights = scene.GetNumLights();
    uint8_t lightIndex = static_cast<uint8_t>(
        glm::min(rng.UniformFloat(), 1.f - glm::epsilon<float>()) * numLights);

    const Light& light = scene.GetLight(lightIndex);
    glm::vec3 f(0.f);

    // Compute direct light
#if BSDF_SAMPLING
    scatteringPdf = 0.f;
    glm::vec2 scatterSample(rng.UniformFloat(), rng.UniformFloat());
    float bsdfSample = rng.UniformFloat();

    f = bsdf.Sample_f(wo, wi, bsdfSample, scatterSample, scatteringPdf, flags,
                      1, eta_outer);

    if (scatteringPdf > 0.f) {
        float flip = wi.z > 0.f ? 1.f : -1.f;
        Ray shadowRay(isect.p + (isect.gn * shadowBias * flip),
                      bsdf.ToWorld(wi));
        Intersection lightIsect;
        Li = light.Li(lightIsect, isect.p, bsdf.ToWorld(wi), &lightingPdf);
        if (!scene.GetBVH().Intersect(shadowRay, lightIsect)) {
            float weight = 1.f;

            if (!(flags & SPECULAR)) {
#if LIGHT_SAMPLING
                weight =
                    (scatteringPdf * scatteringPdf) /
                    (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                if (lightingPdf > 0.f) {
                    L += (f * Li * glm::abs(wi.z) * weight) / scatteringPdf;
                }
            }

            else {
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
    Li =
        light.Sample_Li(lightIsect, isect.p, wiWorld, lightSample, lightingPdf);
    wi = bsdf.ToLocal(wiWorld);

    float flip = wi.z > 0.f ? 1.f : -1.f;
    Ray shadowRay(isect.p + (isect.gn * shadowBias * flip), wiWorld);
    if (!scene.GetBVH().Intersect(shadowRay, lightIsect) && lightingPdf > 0.f) {
        float weight = 1.f;
        scatteringPdf = bsdf.Pdf(wo, wi, 1, eta_outer);
        if (scatteringPdf > 0.f) {
            glm::vec3 f = bsdf.f(wo, wi, 1, eta_outer);
#if BSDF_SAMPLING
            weight =
                (lightingPdf * lightingPdf) /
                (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
            L += (f * Li * glm::abs(wi.z) * weight) / lightingPdf;
        }
    }
#endif

    return L * numLights;
}

void PathIntegrator::UpdateIsectList(IsectInfoList& isectList,
                                     const Intersection& isect,
                                     float eta_sampled) const {
    // Check if intersected mesh already in isectList
    for (uint32_t k = isectList.size(); k-- > 0;) {
        if (isectList[k].meshID == isect.meshID) {
            // If it is, we are exiting it, and it must be
            // removed
            isectList.erase(isectList.begin() + k);
            return;
        }
    }

    // If it isn't, we are entering said mesh, and add
    // it to the list
    isectList.emplace_back(
        IntersectionInfo(isect.meshID, isect.priority, eta_sampled));

    return;
}

glm::vec4 PathIntegrator::Li_alpha(RNG& rng, Ray ray, const Scene& scene,
                                   const RenderParams& params,
                                   MemoryArena& memoryArena) const {
    // List of IDs and priorities of intersections
    IsectInfoList isectList{ArenaAllocator<IntersectionInfo>(&memoryArena)};
    isectList.reserve(params.bounces);

    // Radiance
    glm::vec3 L(0.f, 0.f, 0.f);
    float alpha = 0.f;
    // Currently assuming camera is always sitting inside of a
    // vacuum
    float eta_sampled = 1.f;
    float eta_outer = 1.f;
    Intersection isect;
    // Throughput
    glm::vec3 beta(1.f);
    uint8_t flags = 0;

    float gamma = params.rougheningFactor * params.rougheningFactor;
    float alphaTweak = 1.f;

    for (uint32_t bounce = 0; bounce < params.bounces; ++bounce) {
        // Light intersections
        float lightTMax = isect.tMax;
        bool lightHit = false;
        glm::vec3 Le(0.f);
        for (uint8_t j = 0; j < scene.GetNumLights(); ++j) {
            Intersection lightIsect;
            const Light& light = scene.GetLight(j);
            glm::vec3 Li = light.Li(lightIsect, ray.o, ray.d);
            if (lightIsect.tMax < lightTMax) {
                Le = Li;
                lightTMax = lightIsect.tMax;
                isect.tMax = lightIsect.tMax;
                lightHit = true;
                alpha = 1.f;
            }
        }

        if (scene.Intersect(ray, isect)) {
            BSDF bsdf =
                isect.material->CreateBSDF(isect, alphaTweak, memoryArena);

            if (IsectIsValid(isect, isectList, eta_outer)) {
                if (bounce == 0) alpha = 1.f;

                glm::vec3 wo = bsdf.ToLocal(-ray.d);

                // Estimate direct lighting
                uint8_t directFlags = 0;
                L += EstimateDirect(scene, wo, bsdf, isect, ray, rng,
                                    directFlags, eta_outer) *
                     beta;

                // Spawn new ray
                glm::vec2 scatteringSample(rng.UniformFloat(),
                                           rng.UniformFloat());
                float bsdfSample = rng.UniformFloat();

                float scatteringPdf = 0.f;
                // Sample for new alpha
                float alpha_i;
                glm::vec3 wi;

                glm::vec3 f = bsdf.Sample_f(
                    wo, wi, bsdfSample, scatteringSample, scatteringPdf, flags,
                    0, eta_outer, &alpha_i, &eta_sampled);

                // Sample for new ray
                if (scatteringPdf <= 0.f) break;
                alphaTweak = (1.f - (gamma * alpha_i)) * alphaTweak;
                beta *= (f / scatteringPdf) * glm::abs(wi.z);
                // Transform to world
                float flip = wi.z > 0.f ? 1.f : -1.f;
                ray = Ray(isect.p + (isect.gn * shadowBias * flip),
                          bsdf.ToWorld(wi));
            }

            else {
                ray = Ray(isect.p + ray.d * shadowBias, ray.d);
                flags = TRANSMISSIVE;

                float bsdfSample = rng.UniformFloat();
                eta_sampled = bsdf.Sample_eta(bsdfSample);
            }

            // Add intersection info to list
            if (flags & TRANSMISSIVE) {
                UpdateIsectList(isectList, isect, eta_sampled);
            }

            // Russian roulette
            float q = glm::max((beta.x + beta.y + beta.z) * 0.33333f, 0.f);

            if (bounce > 3) {
                if (q >= rng.UniformFloat()) {
                    beta /= q;
                }

                else
                    break;
            }

            // Reset intersection
            isect = Intersection();
        }

        else if (bounce == 0) {
            if (lightHit == true) {
                L = Le;
            }
            break;
        }
    }
    return glm::vec4(L, alpha);
}