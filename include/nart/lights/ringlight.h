#pragma once

#include "../core/light.h"

class RingLight : public Light
{
public:
    RingLight(float radius, float innerRadius, const glm::vec3& Le, float intensity, const glm::mat4& LightToWorld);

    glm::vec3 Li(Intersection& lightIsect, const glm::vec3& p, const glm::vec3& wi, float* pdf = nullptr) const;

    glm::vec3 Sample_Li(Intersection& lightIsect, const glm::vec3& p, glm::vec3& wi, glm::vec2 sample, float& pdf) const;

    float Pdf(Intersection& lightIsect, const glm::vec3& p, const glm::vec3& wi) const;

private:
    const float innerRadius = 0.f;
    const float radius = 1.f;
};


