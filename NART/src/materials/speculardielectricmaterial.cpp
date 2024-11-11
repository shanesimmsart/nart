#include "speculardielectricmaterial.h"

SpecularDielectricMaterial::SpecularDielectricMaterial(glm::vec3 rho, float eta) : rho(rho), eta(eta)
{}

BSDF SpecularDielectricMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha = 0.f;
    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    /*
    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, glm::max(0.0001f, alpha), alpha_prime));
    }
    */

    // else
    bsdf.AddBxDF(std::make_shared<SpecularDielectricBRDF>(rho, eta));

    return bsdf;
}


