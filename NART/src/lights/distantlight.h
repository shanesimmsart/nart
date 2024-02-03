#pragma once

#include "../core/light.h"

class DistantLight : public Light
{
public:
    DistantLight(glm::vec3 Le, float intensity, glm::mat4 LightToWorld);

    glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi);

    glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf);

    float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi);

    glm::vec3 direction;
};


