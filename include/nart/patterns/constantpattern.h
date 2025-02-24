#pragma once

#include "../core/pattern.h"

class ConstantPattern : public Pattern {
public:
    ConstantPattern(glm::vec3 value);

    glm::vec3 Sample(const glm::vec2& sample, glm::vec2& pdfSample, float& pdf);

    float Pdf(const glm::vec2& sample);

    glm::vec3 GetValue(const Intersection& isect);

private:
    glm::vec3 value = glm::vec3(0.f);
};
