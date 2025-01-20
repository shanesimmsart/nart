#pragma once

#include "../core/material.h"
#include "../core/pattern.h"

class GlassMaterial : public Material {
public:
    GlassMaterial(PatternPtr&& rho_s, PatternPtr&& tau, PatternPtr&& eta,
                  PatternPtr&& alpha);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_sPtn;
    PatternPtr tauPtn;
    PatternPtr etaPtn;
    PatternPtr alphaPtn;
};
