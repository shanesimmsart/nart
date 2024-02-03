#pragma once

#include "../core/material.h"

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    glm::vec3 rho;
    float alpha;
    glm::vec3 R;
    float eta;
};
