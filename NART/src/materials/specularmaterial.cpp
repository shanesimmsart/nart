#include "specularmaterial.h"

SpecularMaterial::SpecularMaterial(glm::vec3 R, float eta) : R(R), eta(eta)
{}

BSDF SpecularMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 1);

    float alpha = glm::min(roughnessOffset, 1.f);

    if (alpha > 0.001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alpha));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}


