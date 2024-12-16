#pragma once

#include "../core/material.h"

class GlossyDielectricMaterial : public Material
{
public:
    GlossyDielectricMaterial(glm::vec3 rho_s, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float alphaTweak);

private:
    float alpha;
    glm::vec3 rho_s;
    float eta;
};


