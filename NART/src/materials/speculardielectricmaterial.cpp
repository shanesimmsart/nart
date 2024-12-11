#include "speculardielectricmaterial.h"

#define SPECULAR 0
#define TS 0

SpecularDielectricMaterial::SpecularDielectricMaterial(glm::vec3 rho, float eta, float alpha) : rho(rho), eta(eta), alpha(alpha)
{}

BSDF SpecularDielectricMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

#if !SPECULAR && !TS
    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_unique<DielectricBRDF>(glm::vec3(1.f), eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_unique<SpecularDielectricBRDF>(rho, eta));
#elif SPECULAR
    bsdf.AddBxDF(std::make_shared<SpecularDielectricBRDF>(rho, eta));
#else
    bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(glm::vec3(1.f), eta, glm::max(0.0001f, alpha), alpha_prime));
#endif

    return bsdf;
}


