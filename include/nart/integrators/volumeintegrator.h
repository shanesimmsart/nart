#pragma once

#include "../core/integrator.h"

class VolumeIntegrator : public Integrator {
public:
    VolumeIntegrator() {}

    glm::vec4 Li_alpha(RNG& rng, Ray ray, const Scene& scene,
                       const RenderParams& params,
                       MemoryArena& memoryArena) const;

private:
    const float shadowBias = 0.001f;
};
