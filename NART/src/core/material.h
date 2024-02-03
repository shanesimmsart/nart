#pragma once

#include "reflection.h"

class Material
{
public:
    virtual BSDF CreateBSDF(glm::vec3 n, float roughnessOffset) = 0;

protected:
    Material() {};
};


