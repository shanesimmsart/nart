#include "../../include/nart/materials/plasticmaterial.h"

PlasticMaterial::PlasticMaterial(PatternPtr&& rho_d, PatternPtr&& rho_s,
                                 PatternPtr&& eta, PatternPtr&& alpha,
                                 PatternPtr&& normalPtn)
    : rho_dPtn(std::move(rho_d)),
      rho_sPtn(std::move(rho_s)),
      etaPtn(std::move(eta)),
      alphaPtn(std::move(alpha)),
      normalPtn(std::move(normalPtn)) {}

BSDF PlasticMaterial::CreateBSDF(const Intersection& isect, float alphaTweak,
                                 MemoryArena& memoryArena) {
    BSDF bsdf(isect, 2);

    if (normalPtn) {
        glm::vec3 n = normalPtn->GetValue(isect);
        n *= 2.f;
        n -= glm::vec3(1.f);
        bsdf.BuildCoordSys(isect, &n);
    }

    else
        bsdf.BuildCoordSys(isect);

    float alpha = alphaPtn->GetValue(isect).x;
    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);
    glm::vec3 rho_d = rho_dPtn->GetValue(isect);
    glm::vec3 rho_s = rho_sPtn->GetValue(isect);
    float eta = etaPtn->GetValue(isect).x;

    BxDF* lambert =
        new (memoryArena.Allocate(sizeof(LambertBRDF))) LambertBRDF(rho_d);
    bsdf.AddBxDF(*lambert);

    if (alpha_prime > 0.001f) {
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
