#pragma once

#include "../core/material.h"

class GlassMaterial : public Material
{
public:
    GlassMaterial(glm::vec3 rho_s, glm::vec3 tau, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float alphaTweak);

private:
    glm::vec3 rho_s;
    glm::vec3 tau;
    float eta;
    float alpha;
};


