#include "diffusematerial.h"

DiffuseMaterial::DiffuseMaterial(glm::vec3 rho) : rho(rho) // your boat...
{}

BSDF DiffuseMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 1);
    bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));
    return bsdf;
}


