#pragma once

#include "../core/material.h"
#include "../core/pattern.h"
#include "../bxdfs/lambertbrdf.h"

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(PatternPtr&& rhoPtn, PatternPtr&& normalPtn = nullptr);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    // Eventually, material inputs will be replaced with patterns, that vary
    // depending on UVs, position, etc.
    PatternPtr rhoPtn;
    PatternPtr normalPtn = nullptr;
};
