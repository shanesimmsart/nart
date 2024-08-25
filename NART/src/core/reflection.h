#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "sampling.h"

#define BECKMANN 0

float Fresnel(float etaI, float etaT, float cosTheta);

enum BSDFFlags
{
    SPECULAR = 1,
    GLOSSY = 2,
    DIFFUSE = 4
};

class BxDF
{
public:
    // BxDF - wo and wi must always be in local coord sys
    virtual glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi) = 0;

    // Sample BRDF, setting wi and pdf
    virtual glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags) = 0;

    virtual float Pdf(const glm::vec3& wo, const glm::vec3& wi) = 0;

    virtual ~BxDF() {}

    BSDFFlags flags;

protected:
    BxDF() {}
};

class BSDF
{
public:
    BSDF(glm::vec3 n, uint8_t numBxDFs);

    glm::vec3 ToLocal(glm::vec3 v)
    {
        return glm::normalize(glm::vec3(glm::dot(v, nt), glm::dot(v, nb), glm::dot(v, n)));
    }

    glm::vec3 ToWorld(glm::vec3 v)
    {
        return glm::normalize(v.x * nt + v.y * nb + v.z * n);
    }

    void AddBxDF(std::shared_ptr<BxDF> bxdf)
    {
        bxdfs.emplace_back(bxdf);
    }

    // Sum BxDFs
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        glm::vec3 f = glm::vec3(0.f);
        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            f += bxdfs[i]->f(wo, wi);
        }
        return f;
    }

    // Sample and sum BRDFs
    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        // Choose a BxDF
        uint8_t bxdfIndex = static_cast<uint8_t>(sample.x * static_cast<float>(numBxDFs));

        // Remap sample to remove bias
        sample.x = glm::fract(sample.x * static_cast<float>(numBxDFs));

        glm::vec3 f = bxdfs[bxdfIndex]->Sample_f(wo, wi, sample, pdf, flags);

        float numEvaluated = 1.f;

        if (!(*flags & SPECULAR))
        {
            for (uint8_t i = 0; i < numBxDFs; ++i)
            {
                if (i != bxdfIndex && !(bxdfs[i]->flags & SPECULAR))
                {
                    float bxdfPdf = bxdfs[i]->Pdf(wo, *wi);
                    if (bxdfPdf > 0.f)
                    {
                        *pdf += bxdfs[i]->Pdf(wo, *wi);
                        f += bxdfs[i]->f(wo, *wi);
                        numEvaluated += 1.f;
                    }
                }
            }
            *pdf /= static_cast<float>(numBxDFs); //numEvaluated;
        }

        return f;
    }

    // Average pdfs
    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        float pdf = 0.f;

        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            pdf += bxdfs[i]->Pdf(wo, wi);
        }

        return pdf / static_cast<float>(numBxDFs);
    }

protected:
    // Coord sys in worldspace
    glm::vec3 nt, nb, n;
    uint8_t numBxDFs = 0;
    std::vector<std::shared_ptr<BxDF>> bxdfs;
};

class LambertBRDF : public BxDF
{
public:
    LambertBRDF(glm::vec3 rho) : rho(rho)
    {
        flags = DIFFUSE;
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        return rho * glm::one_over_pi<float>();
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = DIFFUSE;
        *wi = CosineSampleHemisphere(sample, pdf);
        return f(wo, *wi);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        return wi.z * glm::one_over_pi<float>();
    }

private:
    // Albedo
    glm::vec3 rho;
};

class SpecularBRDF : public BxDF
{
public:
    SpecularBRDF(glm::vec3 R, float eta) : R(R), eta(eta)
    {
        flags = SPECULAR;
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        // Probability of randomly sampling a delta function == 0
        return glm::vec3(0.f);
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = SPECULAR;

        *wi = glm::vec3(-wo.x, -wo.y, wo.z);
        // Delta distribution, so PDF == 1 at this one sample point
        *pdf = 1.f;

        // Mirror-reflection at grazing angles
        if (wi->z == 0.f) return glm::vec3(1.f);

        return (R * Fresnel(1.f, eta, wi->z)) / glm::abs(wi->z);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        return 0.f;
    }

private:
    // Reflectance
    glm::vec3 R;
    // IOR
    float eta;
};

inline glm::vec3 Reflect(const glm::vec3& w1, const glm::vec3& w2)
{
    return (2.f * glm::dot(w1, w2) * w2) - w1;
}

class TorranceSparrowBRDF : public BxDF
{
public:
    TorranceSparrowBRDF(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
    {
        flags = GLOSSY;
    }

    // Masking-shadowing Function (Smith) (Beckmann / GGX)
    float Lambda(const glm::vec3& w)
    {
        float sinTheta = glm::sqrt(1.f - (w.z * w.z));
        float tanTheta = (sinTheta / w.z);
#if BECKMANN
        float a = 1.f / (alpha * tanTheta);
        return (std::erf(a) - 1.f + (glm::exp(-a * a) / (a * glm::sqrt(glm::pi<float>())))) * 0.5f;
#else
        return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) * 0.5f;
#endif
    }

    float G(const glm::vec3& wo, const glm::vec3& wi)
    {
        return 1.f / (1.f + Lambda(wo) + Lambda(wi));
    }

    float G1(const glm::vec3& w)
    {
        return 1.f / (1.f + Lambda(w));
    }

    // Normal Distribution Function (Beckmann / GGX)
    float D(const glm::vec3 wh)
    {
        float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
        float tanTheta = (sinTheta / wh.z);
        float tan2Theta = tanTheta * tanTheta;
        float theta = glm::asin(sinTheta);
        // I think I read somewhere that (x * x) * (x * x) is faster than
        // x * x * x * x
#if BECKMANN
        return (glm::exp(-tan2Theta / (alpha * alpha))) / (glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z)));
#else
        return 1.f / ((glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z))) * (1.f + (tan2Theta / (alpha * alpha))) * (1.f + (tan2Theta / (alpha * alpha))));
#endif
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        if (wo.z < 0.f || wi.z < 0.f) return glm::vec3(0.f);

        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        float g = G(wo, wi);
        float d = D(wh);
        float fr = Fresnel(1.f, eta, glm::dot(wh, wi));

        if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

        return (R * g * d * fr) / (4.f * wo.z * wi.z);
        // return (R  * D(wh)) / (4.f * wo.z * wi.z);
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = SPECULAR;
        if (alpha > 0.001f) *flags = GLOSSY;
        // TODO: PBRT uses 1.62142f as roughness == 1
        // Where does that come from??
        if (alpha >= 1.0f) *flags = DIFFUSE;

#if BECKMANN
        // Sample Beckmann distribution to get wh
        float tanTheta = glm::sqrt(-alpha * alpha * glm::log(1.f - sample.x));
#else
        float tanTheta = glm::sqrt((sample.x * alpha * alpha) / (1.f - sample.x));
#endif
        float phi = 2.f * glm::pi<float>() * sample.y;

        // TODO: There is definitely a better way to do this...
        float theta = glm::atan(tanTheta);
        float x = glm::cos(phi) * glm::sin(theta);
        float y = glm::sin(phi) * glm::sin(theta);
        float z = glm::cos(theta);

        if (z < 0.f)
        {
            *pdf = 0.f;
            return glm::vec3(0.f);
        }

        glm::vec3 wh(x, y, z);

        *wi = Reflect(wo, wh);

        *wi = glm::normalize(*wi);

        *pdf = Pdf(wo, *wi);

        return f(wo, *wi);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        if (wh.z < 0.f) return 0.f;

        // Need to convert 1/dwh to 1/dwi
        float cosThetaH = glm::min(glm::dot(wh, wi), 1.f);
        float d = D(wh);
        return (wh.z * d) / (4.f * cosThetaH);
    }

private:
    // Reflectance
    glm::vec3 R;
    // IOR
    float eta;
    // Roughness
    float alpha;
};

