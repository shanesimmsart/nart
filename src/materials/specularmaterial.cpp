#include "../../include/nart/materials/specularmaterial.h"

SpecularMaterial::SpecularMaterial(PatternPtr&& rho_s, PatternPtr&& eta)
    : rho_sPtn(std::move(rho_s)), etaPtn(std::move(eta)) {}

BSDF SpecularMaterial::CreateBSDF(const Intersection& isect, float alphaTweak,
                                  MemoryArena& memoryArena) {
    BSDF bsdf(isect.sn, 1);

    float alpha = 0.f;
    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);
    glm::vec3 rho_s = rho_sPtn->GetValue(isect);
    float eta = etaPtn->GetValue(isect).x;

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
