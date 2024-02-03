#pragma once

#include "../core/material.h"

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial(glm::vec3 rho);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    // Eventually, material inputs will be replaced with patterns, that vary depending on UVs, position, etc.
    glm::vec3 rho;
};


