#include "plasticmaterial.h"

PlasticMaterial::PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha) : rho(rho), R(R), eta(eta), alpha(alpha)
{}

BSDF PlasticMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 2);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));

    if (alpha_prime > 0.001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}


