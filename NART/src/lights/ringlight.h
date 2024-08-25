#pragma once

#include "../core/light.h"

class RingLight : public Light
{
public:
    RingLight(float radius, float innerRadius, glm::vec3 Le, float intensity, glm::mat4 LightToWorld);

    glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) const;

    glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf) const;

    float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) const;

    float innerRadius = 0.f;
    float radius = 1.f;
};


