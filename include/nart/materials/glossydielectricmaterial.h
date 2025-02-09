#pragma once

#include "../core/material.h"
#include "../core/pattern.h"
#include "../bxdfs/specularbrdf.h"
#include "../bxdfs/torrancesparrowbrdf.h"

class GlossyDielectricMaterial : public Material {
public:
    GlossyDielectricMaterial(PatternPtr&& rho_s, PatternPtr&& eta,
                             PatternPtr&& alpha, PatternPtr&& normalPtn = nullptr);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_sPtn;
    PatternPtr etaPtn;
    PatternPtr alphaPtn;
    PatternPtr normalPtn = nullptr;
};
