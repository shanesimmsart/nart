#pragma once

#include "../core/bxdf.h"

class LambertBRDF : public BxDF
{
public:
    LambertBRDF(const glm::vec3& rho_d);

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime);

    glm::vec3 Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime);

    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime);

private:
    // Diffuse
    const glm::vec3 rho_d;
};


