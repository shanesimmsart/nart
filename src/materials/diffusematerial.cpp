#include "../../include/nart/materials/diffusematerial.h"

DiffuseMaterial::DiffuseMaterial(const glm::vec3& rho)
    : rho(rho)  // your boat...
{}

BSDF DiffuseMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak,
                                 MemoryArena& memoryArena) {
    BSDF bsdf(n, 1);

    BxDF* lambert =
        new (memoryArena.Allocate(sizeof(LambertBRDF))) LambertBRDF(rho);
    bsdf.AddBxDF(*lambert);

    return bsdf;
}
