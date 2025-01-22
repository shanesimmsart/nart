#include "../../include/nart/materials/glossydielectricmaterial.h"

GlossyDielectricMaterial::GlossyDielectricMaterial(PatternPtr&& rho_s,
                                                   PatternPtr&& eta,
                                                   PatternPtr&& alpha)
    : rho_sPtn(std::move(rho_s)),
      etaPtn(std::move(eta)),
      alphaPtn(std::move(alpha)) {}

BSDF GlossyDielectricMaterial::CreateBSDF(const Intersection& isect,
                                          float alphaTweak,
                                          MemoryArena& memoryArena) {
    BSDF bsdf(isect.sn, 1);

    float alpha = alphaPtn->GetValue(isect).x;
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
