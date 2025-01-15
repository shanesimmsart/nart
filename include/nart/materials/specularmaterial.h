#pragma once

#include "../core/material.h"

class SpecularMaterial : public Material {
public:
    SpecularMaterial(const glm::vec3& rho_s, float eta);

    BSDF CreateBSDF(const glm::vec3& n, float alphaTweak,
                    MemoryArena& memoryArena);

private:
    const glm::vec3 rho_s;
    const float eta;
};
