#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>

#include "geometry.h"

class Light
{
public:
    Light(const glm::vec3& Le, float intensity, const glm::mat4& LightToWorld) : Le(Le), intensity(intensity), LightToWorld(LightToWorld) {}

    // Return incident radiance from light, set wi
    virtual glm::vec3 Li(Intersection& lightIsect, const glm::vec3& p, const glm::vec3& wi, float* pdf = nullptr) const = 0;

    // Sample light, return incident radiance, set wi, return pdf
    virtual glm::vec3 Sample_Li(Intersection& lightIsect, const glm::vec3& p, glm::vec3& wi, glm::vec2 sample, float& pdf) const = 0;

    virtual float Pdf(Intersection& lightIsect, const glm::vec3& p, const glm::vec3& wi) const = 0;

    virtual ~Light() {}

    // Does light have a delta distribution?
    bool isDelta = false;

protected:
    // Radiance emitted
    // struct under each child?
    const glm::vec3 Le;
    const float intensity;
    const glm::mat4 LightToWorld;
};

using LightPtr = std::unique_ptr<Light>;


