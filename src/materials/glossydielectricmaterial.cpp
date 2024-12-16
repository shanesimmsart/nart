#include "../../include/nart/materials/glossydielectricmaterial.h"

GlossyDielectricMaterial::GlossyDielectricMaterial(glm::vec3 rho_s, float eta, float alpha) : rho_s(rho_s), eta(eta), alpha(alpha)
{}

BSDF GlossyDielectricMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_unique<TorranceSparrowBRDF>(rho_s, eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_unique<SpecularBRDF>(rho_s, eta));

    return bsdf;
}


