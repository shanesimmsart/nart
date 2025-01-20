#include "../../include/nart/patterns/constantpattern.h"

ConstantPattern::ConstantPattern(glm::vec3 value) : value(value) {}

glm::vec3 ConstantPattern::GetValue(const Intersection& isect) { return value; }