#pragma once

#include "../core/material.h"

class GlassMaterial : public Material {
public:
    GlassMaterial(const glm::vec3& rho_s, const glm::vec3& tau, float eta,
                  float alpha);

    BSDF CreateBSDF(const glm::vec3& n, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    const glm::vec3 rho_s;
    const glm::vec3 tau;
    const float eta;
    float alpha;
};
