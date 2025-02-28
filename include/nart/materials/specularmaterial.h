#pragma once

#include "../bxdfs/specularbrdf.h"
#include "../bxdfs/torrancesparrowbrdf.h"
#include "../core/material.h"
#include "../core/pattern.h"

class SpecularMaterial : public Material {
public:
    SpecularMaterial(PatternPtr&& rho_s, PatternPtr&& eta,
                     PatternPtr&& normalPtn = nullptr);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_sPtn;
    PatternPtr etaPtn;
    PatternPtr normalPtn = nullptr;
};
