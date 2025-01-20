#pragma once

#include "../core/material.h"
#include "../core/pattern.h"

class GlossyDielectricMaterial : public Material {
public:
    GlossyDielectricMaterial(PatternPtr&& rho_s, PatternPtr&& eta,
                             PatternPtr&& alpha);

    BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    PatternPtr rho_sPtn;
    PatternPtr etaPtn;
    PatternPtr alphaPtn;
};
