#pragma once

#include "../core/material.h"
#include "../core/pattern.h"

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(std::unique_ptr<Pattern>&& rho);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    // Eventually, material inputs will be replaced with patterns, that vary
    // depending on UVs, position, etc.
    // const glm::vec3 rho;
    PatternPtr rhoPtn;
};
