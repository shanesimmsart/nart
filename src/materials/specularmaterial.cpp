#include "../../include/nart/materials/specularmaterial.h"

SpecularMaterial::SpecularMaterial(const glm::vec3& rho_s, float eta)
    : rho_s(rho_s), eta(eta) {}

BSDF SpecularMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak,
                                  MemoryArena& memoryArena) {
    BSDF bsdf(n, 1);

    float alpha = 0.f;
    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        BxDF* specularBRDF =
            new (memoryArena.Allocate(sizeof(TorranceSparrowBRDF)))
                TorranceSparrowBRDF(rho_s, eta, glm::max(0.0001f, alpha),
                                    alpha_prime);
        bsdf.AddBxDF(*specularBRDF);
    }

    else {
        BxDF* specularBRDF = new (memoryArena.Allocate(sizeof(SpecularBRDF)))
            SpecularBRDF(rho_s, eta);
        bsdf.AddBxDF(*specularBRDF);
    }

    return bsdf;
}
