#pragma once

#include "geometry.h"

class Pattern {
public:
    virtual glm::vec3 GetValue(const Intersection& isect) = 0;

    virtual ~Pattern() {}

protected:
    Pattern() {}
};

using PatternPtr = std::unique_ptr<Pattern>;


