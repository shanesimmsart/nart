#include "../../include/nart/materials/diffusematerial.h"

DiffuseMaterial::DiffuseMaterial(PatternPtr&& rhoPtn, PatternPtr&& normalPtn)
    : rhoPtn(std::move(rhoPtn)), normalPtn(std::move(normalPtn)) {}

BSDF DiffuseMaterial::CreateBSDF(const Intersection& isect, float alphaTweak,
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

    glm::vec3 rho = rhoPtn->GetValue(isect);

    BxDF* lambert =
        new (memoryArena.Allocate(sizeof(LambertBRDF))) LambertBRDF(rho);
    bsdf.AddBxDF(*lambert);

    return bsdf;
}
