#include "../../include/nart/core/reflection.h"

float Fresnel(float eta_o, float eta_i, float cosTheta)
{
    float cosTheta_o = glm::min(glm::abs(cosTheta), 1.f);

    float sinTheta_o = glm::sqrt(1.f - (cosTheta_o * cosTheta_o));
    float sinTheta_i = (eta_o / eta_i) * sinTheta_o;

    // TIR
    if (sinTheta_i > 1.f) return 1.f;
    float cosTheta_i = glm::sqrt(1.f - (sinTheta_i * sinTheta_i));

    if (glm::abs(cosTheta_o + cosTheta_i) < 0.00001f) return 0.f;

    float fPara = ((eta_i * cosTheta_o) - (eta_o * cosTheta_i)) / ((eta_i * cosTheta_o) + (eta_o * cosTheta_i));
    float fPerp = ((eta_o * cosTheta_o) - (eta_i * cosTheta_i)) / ((eta_o * cosTheta_o) + (eta_i * cosTheta_i));
    return ((fPara * fPara) + (fPerp * fPerp)) * 0.5f;
}

BSDF::BSDF(glm::vec3 n, uint8_t numBxDFs) : n(n), numBxDFs(numBxDFs)
{
    // Build local coord sys from normal
    if (glm::abs(n.x) > glm::abs(n.y))
    {
        n_t = glm::normalize(glm::vec3(n.z, 0.f, -n.x));
    }

    else
    {
        n_t = glm::normalize(glm::vec3(0.f, n.z, -n.y));
    }
    n_b = glm::normalize(glm::cross(n, n_t));

    bxdfs.reserve(numBxDFs);
}

glm::vec3 BSDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    glm::vec3 f = glm::vec3(0.f);
    for (uint8_t i = 0; i < numBxDFs; ++i)
    {
        f += bxdfs[i]->f(wo, wi, use_alpha_prime);
    }
    return f;
}

glm::vec3 BSDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, bool use_alpha_prime, float* alpha_i)
{
    // Choose a BxDF
    uint8_t bxdfIndex = static_cast<uint8_t>(sample.x * static_cast<float>(numBxDFs));

    // rho_semap sample to remove bias
    sample.x = glm::fract(sample.x * static_cast<float>(numBxDFs));

    glm::vec3 f = bxdfs[bxdfIndex]->Sample_f(wo, wi, sample1D, sample, pdf, flags, alpha_i, use_alpha_prime);

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

float BSDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    float pdf = 0.f;

    for (uint8_t i = 0; i < numBxDFs; ++i)
    {
        pdf += bxdfs[i]->Pdf(wo, wi, use_alpha_prime);
    }

    return pdf / static_cast<float>(numBxDFs);
}



LambertBRDF::LambertBRDF(glm::vec3 rho_d) : rho_d(rho_d)
{
    flags = DIFFUSE;
    alpha = 1.f;
    alpha_0 = 1.f;
    alpha_prime = 1.f;
}

glm::vec3 LambertBRDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    // Note: alpha should never exceed 1.0, don't need to adjust
    return rho_d * glm::one_over_pi<float>();
}

glm::vec3 LambertBRDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime)
{
    if (alpha_i != nullptr) *alpha_i = alpha;

    *flags = DIFFUSE;
    *wi = CosineSampleHemisphere(sample, pdf);
    return f(wo, *wi, use_alpha_prime);
}

float LambertBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    return wi.z * glm::one_over_pi<float>();
}



SpecularBRDF::SpecularBRDF(glm::vec3 rho_s, float eta) : rho_s(rho_s), eta(eta)
{
    flags = SPECULAR;
}

glm::vec3 SpecularBRDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    // Probability of randomly sampling a delta function == 0
    return glm::vec3(0.f);
}

glm::vec3 SpecularBRDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime)
{
    if (alpha_i != nullptr) *alpha_i = 0.f;
    *flags = SPECULAR;

    *wi = glm::vec3(-wo.x, -wo.y, wo.z);
    // Delta distribution, so PDF == 1 at this one sample point
    *pdf = 1.f;

    // Mirror-reflection at grazing angles
    if (wi->z == 0.f) return glm::vec3(1.f);

    return (rho_s * Fresnel(1.f, eta, wi->z)) / glm::abs(wi->z);
}

float SpecularBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    return 0.f;
}



SpecularDielectricBRDF::SpecularDielectricBRDF(glm::vec3 rho_s, glm::vec3 tau, float eta) : rho_s(rho_s), tau(tau), eta(eta)
{
    flags = SPECULAR;
}

glm::vec3 SpecularDielectricBRDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    // Probability of randomly sampling a delta function == 0
    return glm::vec3(0.f);
}

glm::vec3 SpecularDielectricBRDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime)
{
    if (alpha_i != nullptr) *alpha_i = 0.f;
    *flags = SPECULAR;

    // Check if inside or outside of medium
    float eta_o = 1.f;
    float eta_i = eta;

    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    float Fr = Fresnel(eta_o, eta_i, glm::abs(wo.z));

    if (sample.x < Fr)
    {
        // return glm::vec3(0.f);

        // Delta distribution, so PDF == 1 at this one sample point
        *pdf = Fr;

        // rho_seflect
        *wi = glm::vec3(-wo.x, -wo.y, wo.z);

        // Mirror-reflection at grazing angles
        if (wi->z == 0.f) return glm::vec3(1.f);

        // Check hemisphere with gn
        return glm::vec3(Fr / glm::abs(wi->z)) * rho_s;
    }

    else
    {
        // return glm::vec3(0.f);
        // Delta distribution, so PDF == 1 at this one sample point
        *pdf = 1.f - Fr;

        // rho_sefract
        float sinTheta_o = glm::sqrt(1.f - (wo.z * wo.z));
        float sinTheta_i = ((eta_o / eta_i) * sinTheta_o);

        // TIR
        if (sinTheta_i >= 1.f)
        {
            // return glm::vec3(0.f);

            *wi = glm::vec3(-wo.x, -wo.y, wo.z);
            return glm::vec3(1.f);
        }

        *flags |= TRANSMISSIVE;

        glm::vec3 a(wo.x, wo.y, 0.f);
        glm::vec3 c = -a * (eta_o / eta_i);
        glm::vec3 d(0.f, 0.f, -glm::sqrt(1.f - (sinTheta_i * sinTheta_i)));
        if (wo.z < 0.f) d = -d;
        *wi = glm::normalize(c + d);

        // Check hemisphere with gn
        glm::vec3 f = glm::vec3(((eta_o / eta_i) * (eta_o / eta_i) * (1.f - Fr)) / glm::abs(wi->z));

        return f;
    }
}

float SpecularDielectricBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    return 0.f;
}



TorranceSparrowBRDF::TorranceSparrowBRDF(glm::vec3 rho_s, float eta, float _alpha_0, float _alpha_prime) : rho_s(rho_s), eta(eta)
{
    alpha = _alpha_0;
    alpha_0 = _alpha_0;
    alpha_prime = _alpha_prime;
    flags = GLOSSY;
}

float TorranceSparrowBRDF::Lambda(const glm::vec3& w)
{
    float sinTheta = glm::sqrt(1.f - (w.z * w.z));
    float tanTheta = (sinTheta / w.z);
    return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) * 0.5f;
}

float TorranceSparrowBRDF::D(const glm::vec3 wh)
{
    float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
    float tanTheta = (sinTheta / wh.z);
    float tan2Theta = tanTheta * tanTheta;
    float theta = glm::asin(sinTheta);
    // I think I read somewhere that (x * x) * (x * x) is faster than
    // x * x * x * x
    return 1.f / ((glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z))) * (1.f + (tan2Theta / (alpha * alpha))) * (1.f + (tan2Theta / (alpha * alpha))));
}

glm::vec3 TorranceSparrowBRDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    if (wo.z < 0.f || wi.z < 0.f) return glm::vec3(0.f);

    glm::vec3 wh = wo + wi;
    wh = glm::normalize(wh);

    float g = G(wo, wi);
    float d = D(wh);
    float fr = Fresnel(1.f, eta, glm::dot(wh, wi));

    if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

    return (rho_s * g * d * fr) / (4.f * wo.z * wi.z);
    // return (rho_s  * D(wh)) / (4.f * wo.z * wi.z);
}

glm::vec3 TorranceSparrowBRDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    if (alpha_i != nullptr) *alpha_i = alpha;
    *flags = SPECULAR;
    if (alpha > 0.001f) *flags = GLOSSY;
    if (alpha >= 1.0f) *flags = DIFFUSE;

    // Transform wo from ellipsoid to hemisphere
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

    *wi = rho_seflect(wo, wh);

    *wi = glm::normalize(*wi);

    *pdf = Pdf(wo, *wi, use_alpha_prime);

    return f(wo, *wi, use_alpha_prime);
}

float TorranceSparrowBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    glm::vec3 wh = wo + wi;
    wh = glm::normalize(wh);

    if (wh.z < 0.f) return 0.f;

    float cosThetaH = glm::min(glm::dot(wo, wh), 1.f);
    float pdf = (D(wh) * glm::min(glm::dot(wo, wh), 1.f) * G1(wo)) / wo.z;
    return glm::max(0.f, pdf / (4.f * cosThetaH));
}



DielectricBRDF::DielectricBRDF(glm::vec3 rho_s, glm::vec3 tau, float eta, float _alpha_0, float _alpha_prime) : rho_s(rho_s), tau(tau), eta(eta)
{
    alpha = _alpha_0;
    alpha_0 = _alpha_0;
    alpha_prime = _alpha_prime;
    flags = GLOSSY;
}

float DielectricBRDF::Lambda(const glm::vec3& w)
{
    float sinTheta = glm::sqrt(1.f - (w.z * w.z));
    float tanTheta = (sinTheta / w.z);
    return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) * 0.5f;
}

float DielectricBRDF::D(const glm::vec3 wh)
{
    if (wh.z == 0.f) return 0.f;
    float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
    float tanTheta = (sinTheta / wh.z);
    float tan2Theta = tanTheta * tanTheta;
    float theta = glm::asin(sinTheta);
    // I think I read somewhere that (x * x) * (x * x) is faster than
    // x * x * x * x
    return 1.f / ((glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z))) * (1.f + (tan2Theta / (alpha * alpha))) * (1.f + (tan2Theta / (alpha * alpha))));
}

glm::vec3 DielectricBRDF::f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    float eta_o = 1.f;
    float eta_i = eta;

    if (wo.z < 0.f) std::swap(eta_o, eta_i);

    if (wo.z * wi.z >= 0.f)
    {
        // return glm::vec3(0.f);

        // Reflect
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);
        if (wh.z < 0.f) wh *= -1.f;

        float g = G(wo, wi);
        float d = D(wh);
        float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));

        if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

        return ((rho_s * g * d * Fr) / (4.f * wo.z * wi.z)) * rho_s;
    }

    else
    {
        // rho_sefract
        glm::vec3 wh = glm::normalize((eta_o * wo) + (eta_i * wi));
        if (wh.z < 0.f) wh *= -1.f;

        float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));
        if (Fr >= 1.f)
        {
            return glm::vec3(0.f);
        }

        float g = G(wo, wi);
        float d = D(wh);
        float wiDotWh = glm::dot(wi, wh);
        float woDotWh = glm::dot(wo, wh);

        float num = g * d * (1.f - Fr) * eta_o * eta_o * glm::abs(wiDotWh) * glm::abs(woDotWh);
        // abs?
        float denom = ((eta_i * wiDotWh) + (eta_o * woDotWh)) * ((eta_i * wiDotWh) + (eta_o * woDotWh)) * glm::abs(wo.z * wi.z);
        return (glm::vec3(num / denom)) * tau;
    }
}

glm::vec3 DielectricBRDF::Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    if (alpha_i != nullptr) *alpha_i = alpha;
    *flags = SPECULAR;
    if (alpha > 0.001f) *flags = GLOSSY;
    if (alpha >= 1.0f) *flags = DIFFUSE;

    // Transform wo from ellipsoid to hemisphere
    glm::vec3 wo_h = glm::vec3(wo.x * alpha, wo.y * alpha, wo.z);
    wo_h = glm::normalize(wo_h);
    if (wo.z < 0.f) wo_h *= -1.f;

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

    // rho_sefract
    float eta_o = 1.f;
    float eta_i = eta;
    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));

    if (sample1D < Fr)
    {
        // return glm::vec3(0.f);

        // rho_seflect
        *wi = rho_seflect(wo, wh);
        *wi = glm::normalize(*wi);

        *pdf = Pdf(wo, *wi, use_alpha_prime) * Fr;

        return f(wo, *wi, use_alpha_prime);
    }

    // return glm::vec3(0.f);

    // rho_sefract
    float cosTheta_o = glm::dot(wo, wh);
    float sinTheta_o = glm::sqrt(1.f - (cosTheta_o * cosTheta_o));
    float sinTheta_i = ((eta_o / eta_i) * sinTheta_o);

    // TIR
    if (sinTheta_i >= 1.f)
    {
        // rho_seflect
        *wi = rho_seflect(wo, wh);
        *wi = glm::normalize(*wi);

        *pdf = Pdf(wo, *wi, use_alpha_prime) * (1.f - Fr);

        return f(wo, *wi, use_alpha_prime);
    }

    *flags |= TRANSMISSIVE;

    glm::vec3 b = wh * cosTheta_o;
    glm::vec3 a = wo - b;
    glm::vec3 c = -a * (eta_o / eta_i);
    glm::vec3 d = -wh * glm::sqrt(1.f - (sinTheta_i * sinTheta_i));
    if (glm::dot(wo, wh) < 0.f) d *= -1.f;
    *wi = glm::normalize(c + d);

    *pdf = Pdf(wo, *wi, use_alpha_prime) * (1.f - Fr);

    return f(wo, *wi, use_alpha_prime);
}

float DielectricBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime)
{
    if (use_alpha_prime) alpha = alpha_prime;
    else alpha = alpha_0;

    if (wo.z * wi.z >= 0.f)
    {
        // rho_seflect
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        if (wh.z < 0.f) wh *= -1.f;

        float cosThetaH = glm::abs(glm::min(glm::dot(wo, wh), 1.f));
        float pdf = (D(wh) * glm::min(glm::dot(wo, wh), 1.f) * G1(wo)) / wo.z;
        return glm::max(0.f, pdf / (4.f * cosThetaH));
    }

    //Refract
    float eta_o = 1.f;
    float eta_i = eta;
    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    glm::vec3 wh = glm::normalize((eta_o * wo) + (eta_i * wi));
    if (wh.z < 0.f) wh *= -1.f;

    float pdf = (D(wh) * glm::min(glm::abs(glm::dot(wo, wh)), 1.f) * G1(wo)) / glm::abs(wo.z);
    float Jdet = (glm::abs(glm::dot(wi, wh)) * eta_i * eta_i) / (((eta_i * glm::dot(wi, wh)) + (eta_o * glm::dot(wo, wh))) * ((eta_i * glm::dot(wi, wh)) + (eta_o * glm::dot(wo, wh))));
    return pdf * Jdet;
}


