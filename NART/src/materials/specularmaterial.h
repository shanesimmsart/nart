#pragma once

#include "../core/material.h"

class SpecularMaterial : public Material
{
public:
    SpecularMaterial(glm::vec3 R, float eta);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    glm::vec3 R;
    float eta;
};


