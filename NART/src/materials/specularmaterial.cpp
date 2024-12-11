#include "specularmaterial.h"

SpecularMaterial::SpecularMaterial(glm::vec3 R, float eta) : R(R), eta(eta)
{}

BSDF SpecularMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha = 0.f;
    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_unique<TorranceSparrowBRDF>(R, eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_unique<SpecularBRDF>(R, eta));

    return bsdf;
}


