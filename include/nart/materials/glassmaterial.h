#pragma once

#include "../bxdfs/dielectricbrdf.h"
#include "../bxdfs/speculardielectricbrdf.h"
#include "../core/material.h"
#include "../core/pattern.h"

class GlassMaterial : public Material {
public:
    GlassMaterial(PatternPtr&& rho_s, PatternPtr&& tau, PatternPtr&& eta,
                  PatternPtr&& alpha, PatternPtr&& normalPtn = nullptr);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_sPtn;
    PatternPtr tauPtn;
    PatternPtr etaPtn;
    PatternPtr alphaPtn;
    PatternPtr normalPtn = nullptr;
};
