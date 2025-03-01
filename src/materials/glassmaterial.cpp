#include "../../include/nart/materials/glassmaterial.h"

GlassMaterial::GlassMaterial(PatternPtr&& rho_s, PatternPtr&& tau,
                             PatternPtr&& eta, PatternPtr&& alpha, PatternPtr&&)
    : rho_sPtn(std::move(rho_s)),
      tauPtn(std::move(tau)),
      etaPtn(std::move(eta)),
      alphaPtn(std::move(alpha)),
      normalPtn(std::move(normalPtn)) {}

BSDF GlassMaterial::CreateBSDF(const Intersection& isect, float alphaTweak,
                               MemoryArena& memoryArena) {
    BSDF bsdf(isect, 1);

    if (normalPtn) {
        glm::vec3 n = normalPtn->GetValue(isect);
        n *= 2.f;
        n -= glm::vec3(1.f);
        bsdf.BuildCoordSys(isect, &n);
    }

    else
        bsdf.BuildCoordSys(isect);

    float alpha = alphaPtn->GetValue(isect).x;
    float alpha_prime =
        1.f - ((1.f - alphaPtn->GetValue(isect).x) * alphaTweak);
    glm::vec3 rho_s = rho_sPtn->GetValue(isect);
    glm::vec3 tau = tauPtn->GetValue(isect);
    float eta = etaPtn->GetValue(isect).x;

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
