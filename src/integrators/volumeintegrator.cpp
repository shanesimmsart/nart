#include "../../include/nart/integrators/volumeintegrator.h"

glm::vec4 VolumeIntegrator::Li_alpha(RNG& rng, Ray ray, const Scene& scene,
                                     const RenderParams& params,
                                     MemoryArena& memoryArena) const {
    // Radiance
    glm::vec3 L(0.f, 0.f, 0.f);
    // Throughput
    glm::vec3 beta(1.f);
    uint32_t bounce = 0;
    float alpha = 1.f;

    while (true) {
        bool scattered = false, terminated = false;
        // Sample using T_maj's distribution
        // callback function chooses absorb / scatter / null
        // callback() must have the following signature:
        // bool callback(glm::vec3 p, MediumProperties mp, float T_maj, float
        // sigma_maj);
        float u = rng.UniformFloat();
        float uMode = rng.UniformFloat();
        SampleT_maj(u, ray, rng,
                    [&](glm::vec3 p, MediumProperties mp, float T_maj,
                        float sigma_maj) {
                        // callback function
                        float pAbsorb = mp.sigma_a / sigma_maj;
                        float pScatter = mp.sigma_s / sigma_maj;
                        // float pNull = glm::max(0.f, 1.f - (pAbsorb +
                        // pScatter));
                        if (uMode < pAbsorb) {
                            // absorb
                            terminated = true;
                            L += mp.Le * beta;
                            return false;
                        } else if (uMode < pAbsorb + pScatter) {
                            // scatter
                            if (bounce++ > params.bounces) {
                                terminated = true;
                                return false;
                            }
                            ray.o = p;

                            // Isotropic media only for now
                            glm::vec2 sample2D(rng.UniformFloat(),
                                               rng.UniformFloat());
                            float pdf = 0.f;
                            mp.pf.SamplePhaseFunction(sample2D, ray.d, pdf);
                            // Current phase function = it's own PDF
                            // beta *= mp.pf / pdf;

                            scattered = true;
                            return false;
                        } else {
                            // null scatter
                            uMode = rng.UniformFloat();
                            return true;
                        }
                    });
        if (terminated) {
            break;
        }
        if (scattered) {
            continue;
        }

        // Intersect lights
        float lightTMax = Infinity;
        bool lightHit = false;
        glm::vec3 Le(0.f);
        for (uint8_t j = 0; j < scene.GetNumLights(); ++j) {
            Intersection lightIsect;
            const Light& light = scene.GetLight(j);
            glm::vec3 Li = light.Li(lightIsect, ray.o, ray.d);
            if (lightIsect.tMax < lightTMax) {
                Le = Li;
                lightTMax = lightIsect.tMax;
            }
        }
        L += Le * beta;
        break;
    }

    return glm::vec4(L, alpha);
}
