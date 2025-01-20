#pragma once

#include "../core/material.h"
#include "../core/pattern.h"

class PlasticMaterial : public Material {
public:
    PlasticMaterial(PatternPtr&& rho_d, PatternPtr&& rho_s, PatternPtr&& eta,
                    PatternPtr&& alpha);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_dPtn;
    PatternPtr rho_sPtn;
    PatternPtr etaPtn;
    PatternPtr alphaPtn;
};
