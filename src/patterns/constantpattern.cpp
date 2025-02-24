#include "../../include/nart/patterns/constantpattern.h"

ConstantPattern::ConstantPattern(glm::vec3 value) : value(value) {}

glm::vec3 ConstantPattern::Sample(const glm::vec2& sample, glm::vec2& pdfSample,
                                  float& pdf) {
    pdfSample = sample;
    pdf = 1.f;
    return value;
}

float ConstantPattern::Pdf(const glm::vec2& sample) { return 1.f; }

glm::vec3 ConstantPattern::GetValue(const Intersection& isect) { return value; }