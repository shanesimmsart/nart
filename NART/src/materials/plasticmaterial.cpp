#include "plasticmaterial.h"

PlasticMaterial::PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha) : rho(rho), R(R), eta(eta), alpha(alpha)
{}

BSDF PlasticMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 2);

    bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));

    float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

    if (alphaAdjusted > 0.001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}


