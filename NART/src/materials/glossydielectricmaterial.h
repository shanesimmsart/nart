#pragma once

#include "../core/material.h"

class GlossyDielectricMaterial : public Material
{
public:
    GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    float alpha;
    glm::vec3 R;
    float eta;
};


