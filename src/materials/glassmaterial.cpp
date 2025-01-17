#include "../../include/nart/materials/glassmaterial.h"

GlassMaterial::GlassMaterial(const glm::vec3& rho_s, const glm::vec3& tau,
                             float eta, float alpha)
    : rho_s(rho_s), tau(tau), eta(eta), alpha(alpha) {}

BSDF GlassMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak,
                               MemoryArena& memoryArena) {
    BSDF bsdf(n, 1);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        BxDF* lambert =
            new (memoryArena.Allocate(sizeof(DielectricBRDF))) DielectricBRDF(
                rho_s, tau, eta, glm::max(0.0001f, alpha), alpha_prime);
        bsdf.AddBxDF(*lambert);
    }

    else {
        BxDF* specularBRDF =
            new (memoryArena.Allocate(sizeof(SpecularDielectricBRDF)))
                SpecularDielectricBRDF(rho_s, tau, eta);
        bsdf.AddBxDF(*specularBRDF);
    }

    return bsdf;
}
