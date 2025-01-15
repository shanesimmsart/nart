#pragma once

#include "../core/material.h"

class GlossyDielectricMaterial : public Material {
public:
    GlossyDielectricMaterial(const glm::vec3& rho_s, float eta, float alpha);

    BSDF CreateBSDF(const glm::vec3& n, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    const glm::vec3 rho_s;
    const float eta;
    float alpha;
};
