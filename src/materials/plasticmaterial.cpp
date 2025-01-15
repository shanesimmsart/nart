#include "../../include/nart/materials/plasticmaterial.h"

PlasticMaterial::PlasticMaterial(const glm::vec3& rho_d, const glm::vec3& rho_s, float eta, float alpha) : rho_d(rho_d), rho_s(rho_s), eta(eta), alpha(alpha)
{}

BSDF PlasticMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak, MemoryArena& memoryArena)
{
    BSDF bsdf(n, 2);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    BxDF* lambert = new(memoryArena.Allocate(sizeof(LambertBRDF))) LambertBRDF(rho_d);
    bsdf.AddBxDF(*lambert);

    if (alpha_prime > 0.001f) {
        BxDF* specularBRDF = new(memoryArena.Allocate(sizeof(TorranceSparrowBRDF))) TorranceSparrowBRDF(rho_s, eta, glm::max(0.0001f, alpha), alpha_prime);
        bsdf.AddBxDF(*specularBRDF);
    }

    else
    {
        BxDF* specularBRDF = new(memoryArena.Allocate(sizeof(SpecularBRDF))) SpecularBRDF(rho_s, eta);
        bsdf.AddBxDF(*specularBRDF);
    }

    return bsdf;
}


