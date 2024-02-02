#include "materials.h"

DiffuseMaterial::DiffuseMaterial(glm::vec3 rho) : rho(rho) // your boat...
{}

BSDF DiffuseMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 1);
    bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));
    return bsdf;
}



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



GlossyDielectricMaterial::GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
{}

BSDF GlossyDielectricMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 1);

    float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

    if (alphaAdjusted > 0.001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}



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


