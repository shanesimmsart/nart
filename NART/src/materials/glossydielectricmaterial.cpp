#include "glossydielectricmaterial.h"

GlossyDielectricMaterial::GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
{}

BSDF GlossyDielectricMaterial::CreateBSDF(glm::vec3 n, float roughnessOffset)
{
    BSDF bsdf(n, 1);

    float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

    if (alphaAdjusted > 0.001f) {
        bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
    }

    else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

    return bsdf;
}


