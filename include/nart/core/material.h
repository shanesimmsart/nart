#pragma once

#include "../bxdfs/dielectricbrdf.h"
#include "../bxdfs/lambertbrdf.h"
#include "../bxdfs/specularbrdf.h"
#include "../bxdfs/speculardielectricbrdf.h"
#include "../bxdfs/torrancesparrowbrdf.h"

class Material
{
public:
    // Creates BSDF and adds BxDFs
    virtual BSDF CreateBSDF(const glm::vec3& n, float alphaTweak) = 0;

    virtual ~Material() {}

protected:
    Material() {};
};

using MaterialPtr = std::unique_ptr<Material>;


