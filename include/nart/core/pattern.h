#pragma once

#include "bxdf.h"

class Pattern {
public:
    virtual glm::vec3 GetValue(const Intersection& isect) = 0;

    virtual glm::vec3 Sample(const glm::vec2& sample, glm::vec2& pdfSample,
                             float& pdf) = 0;

    virtual float Pdf(const glm::vec2& sample) = 0;

    virtual ~Pattern() {}

protected:
    Pattern() {}
};

using PatternPtr = std::unique_ptr<Pattern>;
