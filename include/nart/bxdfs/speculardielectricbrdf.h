#pragma once

#include "../core/bxdf.h"

// Perfectly specular transmissive BRDF with delta distribution
class SpecularDielectricBRDF : public BxDF {
public:
    SpecularDielectricBRDF(const glm::vec3& rho_s, const glm::vec3& tau,
                           float eta);

    // Note: We should never be in a situation where alpha_prime > 0 and we are
    // using a specular BRDF
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime);

    glm::vec3 Sample_f(const glm::vec3& wo, glm::vec3& wi, float sample1D,
                       glm::vec2 sample, float& pdf, uint8_t& flags,
                       float* alpha_i, bool use_alpha_prime);

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime);

private:
    // Reflectance
    const glm::vec3 rho_s;
    // Transmittance
    const glm::vec3 tau;
    // IOR
    float eta;
};