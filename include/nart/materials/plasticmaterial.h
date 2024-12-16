#pragma once

#include "../core/material.h"

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(glm::vec3 rho_d, glm::vec3 rho_s, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float alphaTweak);

private:
    glm::vec3 rho_d;
    float alpha;
    glm::vec3 rho_s;
    float eta;
};
