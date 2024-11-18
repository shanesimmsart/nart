#pragma once

#include "../core/material.h"

class SpecularDielectricMaterial : public Material
{
public:
    SpecularDielectricMaterial(glm::vec3 rho, float eta);

    BSDF CreateBSDF(glm::vec3 n, float alphaTweak);

private:
    glm::vec3 rho;
    float eta;
};


