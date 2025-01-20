#include "../../include/nart/materials/diffusematerial.h"

DiffuseMaterial::DiffuseMaterial(PatternPtr&& rho)
    : rhoPtn(std::move(rho))
{}

BSDF DiffuseMaterial::CreateBSDF(const Intersection& isect, float alphaTweak,
                                 MemoryArena& memoryArena) {
    BSDF bsdf(isect.sn, 1);

    glm::vec3 rho = rhoPtn->GetValue(isect);

    BxDF* lambert =
        new (memoryArena.Allocate(sizeof(LambertBRDF))) LambertBRDF(rho);
    bsdf.AddBxDF(*lambert);

    return bsdf;
}
