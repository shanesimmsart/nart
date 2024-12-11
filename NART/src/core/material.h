#pragma once

#include "reflection.h"

class Material
{
public:
    // Creates BSDF and adds BxDFs
    virtual BSDF CreateBSDF(glm::vec3 n, float alphaTweak) = 0;

    virtual ~Material() {}

protected:
    Material() {};
};


