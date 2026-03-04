#pragma once

#include "../core/integrator.h"

class PathIntegrator : public Integrator {
public:
    PathIntegrator() {}

    struct IntersectionInfo {
        IntersectionInfo(uint32_t meshID, uint8_t priority, float eta)
            : meshID(meshID), priority(priority), eta(eta) {}

        uint32_t meshID = -1;
        uint8_t priority = 0;
        float eta = 1.f;
    };

    using IsectInfoList =
        std::vector<IntersectionInfo, ArenaAllocator<IntersectionInfo>>;

    // Check if we should ignore the intersection and update IOR of current
    // medium
    bool IsectIsValid(const Intersection& isect, const IsectInfoList& isectList,
                      float& eta_outer) const;

    void UpdateIsectList(IsectInfoList& isectList, const Intersection& isect,
                         float eta_sampled) const;

    glm::vec3 EstimateDirect(const Scene& scene, const glm::vec3 wo, BSDF& bsdf,
                             const Intersection& isect, const Ray& ray,
                             RNG& rng, uint8_t& flags, float eta_outer) const;

    glm::vec4 Li_alpha(RNG& rng, Ray ray, const Scene& scene,
                       const RenderParams& params,
                       MemoryArena& memoryArena) const;

private:
    const float shadowBias = 0.001f;
};