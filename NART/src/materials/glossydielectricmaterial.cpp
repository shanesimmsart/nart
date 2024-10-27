#include "glossydielectricmaterial.h"

GlossyDielectricMaterial::GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
{}

BSDF GlossyDielectricMaterial::CreateBSDF(glm::vec3 n, float alphaTweak)
{
    BSDF bsdf(n, 1);

    float alpha_prime = 1.f - ((1.f - alpha) * alphaTweak);

    if (alpha_prime > 0.0001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, glm::max(0.0001f, alpha), alpha_prime));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}


