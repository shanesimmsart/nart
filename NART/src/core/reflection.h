#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "sampling.h"

#define VNDF 1

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
    virtual glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime) = 0;

    // Sample BRDF, setting wi and pdf
    virtual glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags, float* alpha_i, bool use_alpha_prime) = 0;

    virtual float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime) = 0;

    virtual ~BxDF() {}

    BSDFFlags flags;

protected:
    BxDF() {}

    float alpha = 0.f;
    float alpha0 = 0.f;
    float alpha_prime = 0.f;
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
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        glm::vec3 f = glm::vec3(0.f);
        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            f += bxdfs[i]->f(wo, wi, use_alpha_prime);
        }
        return f;
    }

    // Sample and sum BRDFs
    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags, float* alpha_i, bool use_alpha_prime)
    {
        // Choose a BxDF
        uint8_t bxdfIndex = static_cast<uint8_t>(sample.x * static_cast<float>(numBxDFs));

        // Remap sample to remove bias
        sample.x = glm::fract(sample.x * static_cast<float>(numBxDFs));

        glm::vec3 f = bxdfs[bxdfIndex]->Sample_f(wo, wi, sample, pdf, flags, alpha_i, use_alpha_prime);

        float numEvaluated = 1.f;

        if (!(*flags & SPECULAR))
        {
            for (uint8_t i = 0; i < numBxDFs; ++i)
            {
                if (i != bxdfIndex && !(bxdfs[i]->flags & SPECULAR))
                {
                    float bxdfPdf = bxdfs[i]->Pdf(wo, *wi, use_alpha_prime);
                    if (bxdfPdf > 0.f)
                    {
                        *pdf += bxdfs[i]->Pdf(wo, *wi, use_alpha_prime);
                        f += bxdfs[i]->f(wo, *wi, use_alpha_prime);
                        numEvaluated += 1.f;
                    }
                }
            }
            *pdf /= static_cast<float>(numBxDFs); //numEvaluated;
        }

        return f;
    }

    // Average pdfs
    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        float pdf = 0.f;

        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            pdf += bxdfs[i]->Pdf(wo, wi, use_alpha_prime);
        }

        return pdf / static_cast<float>(numBxDFs);
    }

protected:
    // Coord sys in worldspace
    glm::vec3 nt, nb, n;
    uint8_t numBxDFs = 0;
    // TODO: Look into raw const pointers?
    std::vector<std::shared_ptr<BxDF>> bxdfs;
};

class LambertBRDF : public BxDF
{
public:
    LambertBRDF(glm::vec3 rho) : rho(rho)
    {
        flags = DIFFUSE;
        alpha = 1.f;
        alpha0 = 1.f;
        alpha_prime = 1.f;
    }

    // Note: alpha should never exceed 1.0, don't need to adjust
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        return rho * glm::one_over_pi<float>();
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags, float* alpha_i, bool use_alpha_prime)
    {
        *alpha_i = alpha;
        *flags = DIFFUSE;
        *wi = CosineSampleHemisphere(sample, pdf);
        return f(wo, *wi, use_alpha_prime);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
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

    // Note: We should never be in a situation where alpha_prime > 0 and we are using a specular BRDF
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        // Probability of randomly sampling a delta function == 0
        return glm::vec3(0.f);
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags, float* alpha_i, bool use_alpha_prime)
    {
        *alpha_i = 0.f;
        *flags = SPECULAR;

        *wi = glm::vec3(-wo.x, -wo.y, wo.z);
        // Delta distribution, so PDF == 1 at this one sample point
        *pdf = 1.f;

        // Mirror-reflection at grazing angles
        if (wi->z == 0.f) return glm::vec3(1.f);

        return (R * Fresnel(1.f, eta, wi->z)) / glm::abs(wi->z);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
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
    TorranceSparrowBRDF(glm::vec3 R, float eta, float _alpha0, float _alpha_prime) : R(R), eta(eta)
    {
        alpha = _alpha0;
        alpha0 = _alpha0;
        alpha_prime = _alpha_prime;
        flags = GLOSSY;
    }

    // Masking-shadowing Function (Smith) (Beckmann / GGX)
    float Lambda(const glm::vec3& w)
    {
        float sinTheta = glm::sqrt(1.f - (w.z * w.z));
        float tanTheta = (sinTheta / w.z);
        return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) * 0.5f;
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
        return 1.f / ((glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z))) * (1.f + (tan2Theta / (alpha * alpha))) * (1.f + (tan2Theta / (alpha * alpha))));
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        if (use_alpha_prime) alpha = alpha_prime;
        else alpha = alpha0;

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

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags, float* alpha_i, bool use_alpha_prime)
    {
        if (use_alpha_prime) alpha = alpha_prime;
        else alpha = alpha0;

        *alpha_i = alpha;
        *flags = SPECULAR;
        if (alpha > 0.001f) *flags = GLOSSY;
        // TODO: PBRT uses 1.62142f as roughness == 1
        // Where does that come from??
        if (alpha >= 1.0f) *flags = DIFFUSE;

#if VNDF
        // Transform wi from ellipsoid to hemisphere
        glm::vec3 wo_h = glm::vec3(wo.x * alpha, wo.y * alpha, wo.z);
        wo_h = glm::normalize(wo_h);

        // Sample projection of hemisphere in wo
        // Build coord sys
        glm::vec3 T1 = glm::vec3(wo_h.y, -wo_h.x, 0.f); // cross(wo, z)
        T1 = glm::normalize(T1);
        glm::vec3 T2 = glm::cross(T1, wo_h);
        T2 = glm::normalize(T2);

        // Sample disk
        glm::vec2 visibleHSample = UniformSampleDisk(sample);
        // Convert to visible hemisphere coords
        float s = (1.f + wo_h.z) * 0.5f;
        visibleHSample.y = (s * visibleHSample.y) + ((1.f - s) * glm::sqrt(1.f - (visibleHSample.x * visibleHSample.x)));
        // Project onto hemisphere
        glm::vec3 wh = glm::vec3(glm::sqrt(1.f - (visibleHSample.x * visibleHSample.x) - (visibleHSample.y * visibleHSample.y)), visibleHSample.x, visibleHSample.y);
        // Convert to local coords
        wh = wh.x * wo_h + wh.y * T1 + wh.z * T2;

        // Transform wh to ellipsoid (Inverse transpose)
        wh = glm::vec3(wh.x * alpha, wh.y * alpha, wh.z);
        wh = glm::normalize(wh);
#else
        float tanTheta = glm::sqrt((sample.x * alpha * alpha) / (1.f - sample.x));
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

#endif
        *wi = Reflect(wo, wh);

        *wi = glm::normalize(*wi);

        *pdf = Pdf(wo, *wi, use_alpha_prime);

        return f(wo, *wi, use_alpha_prime);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
    {
        if (use_alpha_prime) alpha = alpha_prime;
        else alpha = alpha0;

        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        if (wh.z < 0.f) return 0.f;

#if VNDF
        float cosThetaH = glm::min(glm::dot(wo, wh), 1.f);
        float pdf = (D(wh) * glm::min(glm::dot(wo, wh), 1.f) * G1(wo)) / wo.z;
        return glm::max(0.f, pdf / (4.f * cosThetaH));
#else
        // Need to convert 1/dwh to 1/dwi
        float cosThetaH = glm::min(glm::dot(wh, wi), 1.f);
        float d = D(wh);
        return (wh.z * d) / (4.f * cosThetaH);
#endif
    }

private:
    // Reflectance
    glm::vec3 R;
    // IOR
    float eta;
};

