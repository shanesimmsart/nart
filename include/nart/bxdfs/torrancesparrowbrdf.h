#pragma once

#include "../core/bxdf.h"

// Microfacet reflective BRDF with Trowbridge-Reitz distribution
class TorranceSparrowBRDF : public BxDF {
public:
    TorranceSparrowBRDF(const glm::vec3& rho_s, float eta, float _alpha_0,
                        float _alpha_prime);

    // Masking-shadowing Function (Smith)
    float Lambda(const glm::vec3& w) const;

    inline float G(const glm::vec3& wo, const glm::vec3& wi) const {
        return 1.f / (1.f + Lambda(wo) + Lambda(wi));
    }

    inline float G1(const glm::vec3& w) const {
        return 1.f / (1.f + Lambda(w));
    }

    // Normal Distribution Function (Trowbridge-Reitz)
    float D(const glm::vec3 wh) const;

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime,
                float eta_outer);

    glm::vec3 Sample_f(const glm::vec3& wo, glm::vec3& wi, float sample1D,
                       glm::vec2 sample, float& pdf, uint8_t& flags,
                       float* alpha_i, bool use_alpha_prime, float eta_outer);

    float Get_eta() const;

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime,
              float eta_outer);

private:
    // alpha used when calculating BRDF
    float alpha = 0.f;
    // unadjusted alpha
    float alpha_0 = 0.f;
    // alpha adjusted via roughening over paths
    float alpha_prime = 0.f;

    // Reflectance
    const glm::vec3 rho_s;
    // IOR
    const float eta;
};
