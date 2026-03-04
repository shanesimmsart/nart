#pragma once

#include "geometry.h"
#include "light.h"
#include "media.h"
#include "scene.h"

class Integrator {
public:
    virtual ~Integrator() {}

    virtual glm::vec4 Li_alpha(RNG& rng, Ray ray, const Scene& scene,
                               const RenderParams& params,
                               MemoryArena& memoryArena) const = 0;

protected:
    Integrator(){};
};

using IntegratorPtr = std::unique_ptr<Integrator>;
