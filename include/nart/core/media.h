#pragma once

#include <array>
#include <utility>

#include "geometry.h"
#include "rng.h"
#include "sampling.h"

struct DensityGrid {
    DensityGrid(uint8_t resolutionX, uint8_t resolutionY, uint8_t resolutionZ,
                std::vector<float> density)
        : resolutionX(resolutionX),
          resolutionY(resolutionY),
          resolutionZ(resolutionZ),
          density(std::move(density)) {}

    float LookUp(uint8_t x, uint8_t y, uint8_t z) const;

    // p is within [0, 1)^3
    float LookUp(glm::vec3 p) const;

    uint8_t resolutionX = 0, resolutionY = 0, resolutionZ = 0;
    std::vector<float> density;
};

using DensityGridPtr = std::unique_ptr<DensityGrid>;

struct MajorantGrid {
    MajorantGrid(const DensityGrid& densityGrid, float sigma_maj);

    float LookUp(uint8_t x, uint8_t y, uint8_t z) const;

    // MajorantGrid is always 16^3
    static constexpr uint8_t width = 1;
    static constexpr uint16_t size = width * width * width;
    static constexpr float invWidth = 1.f / static_cast<float>(width);
    std::array<float, size> majorants;
    float sigma_maj = 0.f;
};

using MajorantGridPtr = std::unique_ptr<MajorantGrid>;

struct RayMajorantSegment {
    RayMajorantSegment() {}

    RayMajorantSegment(float sigma_maj, float tMin, float tMax)
        : sigma_maj(sigma_maj), tMin(tMin), tMax(tMax) {}

    float sigma_maj = 0.f;
    float tMin = Infinity;
    float tMax = 0.f;
};

class RayMajorantIterator {
public:
    RayMajorantIterator() {}

    RayMajorantIterator(const Ray& ray, const glm::vec3& boundsMin,
                        const glm::vec3& boundsMax, MajorantGrid* majorantGrid,
                        float tMin, float tMax);

    bool Next(RayMajorantSegment& seg);

private:
    // majorantGrid has life time of Medium
    MajorantGrid* majorantGrid = nullptr;
    float tMin = Infinity;
    float tCurrent = Infinity;
    float tMax = 0.f;
    uint32_t currentIndex, finalIndex;
    glm::vec3 nextCrossing;
    glm::vec3 crossDistance;
    glm::ivec3 stepAxis;
    bool done = false;
    static constexpr float rayBias = 0.0001f;
};

class PhaseFunction {
public:
    PhaseFunction() {}

    float SamplePhaseFunction(glm::vec2 sample2D, glm::vec3& wi, float& pdf);
};

struct MediumProperties {
    MediumProperties() {}

    MediumProperties(PhaseFunction pf, float sigma_a, float sigma_s,
                     glm::vec3 Le)
        : pf(pf), sigma_a(sigma_a), sigma_s(sigma_s), Le(Le) {}

    PhaseFunction pf;
    float sigma_a = 0.f, sigma_s = 0.f;
    glm::vec3 Le = glm::vec3(0.f);
};

class Medium {
public:
    Medium(glm::vec3 boundsMin, glm::vec3 boundsMax, float sigma_a,
           float sigma_s, glm::vec3 Le, DensityGrid densityGrid)
        : boundsMin(boundsMin),
          boundsMax(boundsMax),
          sigma_a(sigma_a),
          sigma_s(sigma_s),
          Le(Le),
          densityGrid(densityGrid),
          majorantGrid(MajorantGrid(densityGrid, sigma_a + sigma_s)) {}

    bool SampleMedium(glm::vec3(p), MediumProperties& mp);

    bool SampleRay(const Ray& ray, RayMajorantIterator& iter);

private:
    float sigma_a, sigma_s;
    glm::vec3 Le;
    DensityGrid densityGrid;
    MajorantGrid majorantGrid;
    glm::vec3 boundsMin, boundsMax;
};

using MediumPtr = std::unique_ptr<Medium>;

// callback() must have the following signature:
// bool callback(glm::vec3 p, MediumProperties mp, float T_maj, float
// sigma_maj);

template <typename F>
void SampleT_maj(float u, Ray ray, RNG& rng, F callback) {
    RayMajorantIterator iter;
    if (!ray.medium) return;
    if (!ray.medium->SampleRay(ray, iter)) return;

    bool done = false;
    float T_maj = 1.f;

    while (!done) {
        // Attempt to sample segment
        RayMajorantSegment seg;
        if (!iter.Next(seg)) {
            done = true;
            break;
        }
        float tMin = seg.tMin;

        while (true) {
            // Attempt to sample within segment
            float t = tMin +
                      SampleExponentialDecay(rng.UniformFloat(), seg.sigma_maj);
            // u = rng.UniformFloat();
            if (t < seg.tMax) {
                // Check for a real scattering event
                // True means null scattering, keep going
                glm::vec3 p = ray.o + (ray.d * t);
                MediumProperties mp;
                if (!ray.medium->SampleMedium(p, mp)) {
                    // We have left the medium
                    done = true;
                    break;
                }
                if (!callback(p, mp, T_maj, seg.sigma_maj)) {
                    // Real scattering event
                    float dt = t - tMin;
                    T_maj *= glm::exp(-seg.sigma_maj * dt);
                    done = true;
                    break;
                }

                // If we keep going, start sampling from current point
                T_maj = 1.f;
                tMin = t;
            }

            else {
                float dt = seg.tMax - seg.tMin;
                T_maj *= glm::exp(-seg.sigma_maj * dt);
                break;
            }
        }
    }
}
