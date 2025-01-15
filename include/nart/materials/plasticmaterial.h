#pragma once

#include "../core/material.h"

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(const glm::vec3& rho_d, const glm::vec3& rho_s, float eta, float alpha);

    BSDF CreateBSDF(const glm::vec3& n, float alphaTweak, MemoryArena& memoryArena);

private:
    const glm::vec3 rho_d;
    const glm::vec3 rho_s;
    const float eta;
    float alpha;
};
