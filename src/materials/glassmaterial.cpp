#include "../../include/nart/materials/glassmaterial.h"

GlassMaterial::GlassMaterial(const glm::vec3& rho_s, const glm::vec3& tau, float eta, float alpha) : rho_s(rho_s), tau(tau), eta(eta), alpha(alpha)
{}

BSDF GlassMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_unique<DielectricBRDF>(rho_s, tau, eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_unique<SpecularDielectricBRDF>(rho_s, tau, eta));

    return bsdf;
}


