#pragma once

#include "reflection.h"

class Material
{
public:
    virtual BSDF CreateBSDF(glm::vec3 n, float roughnessOffset) = 0;

protected:
    Material() {};
};

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial(glm::vec3 rho) : rho(rho) // your boat...
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);
        bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));
        return bsdf;
    }

private:
    // Eventually, material inputs will be replaced with patterns, that vary depending on UVs, position, etc.
    glm::vec3 rho;
};

class SpecularMaterial : public Material
{
public:
    SpecularMaterial(glm::vec3 R, float eta) : R(R), eta(eta)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);

        float alpha = glm::min(roughnessOffset, 1.f);

        if (alpha > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alpha));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    glm::vec3 R;
    float eta;
};

class GlossyDielectricMaterial : public Material
{
public:
    GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);

        float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

        if (alphaAdjusted > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    float alpha;
    glm::vec3 R;
    float eta;
};

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha) : rho(rho), R(R), eta(eta), alpha(alpha)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 2);

        bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));

        float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

        if (alphaAdjusted > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    glm::vec3 rho;
    float alpha;
    glm::vec3 R;
    float eta;
};


