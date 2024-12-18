#include "../../include/nart/materials/diffusematerial.h"

DiffuseMaterial::DiffuseMaterial(const glm::vec3& rho) : rho(rho) // your boat...
{}

BSDF DiffuseMaterial::CreateBSDF(const glm::vec3& n, float alphaTweak)
{
    BSDF bsdf(n, 1);
    bsdf.AddBxDF(std::make_unique<LambertBRDF>(rho));
    return bsdf;
}


