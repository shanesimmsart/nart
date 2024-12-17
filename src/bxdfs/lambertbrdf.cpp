#include "../../include/nart/bxdfs/lambertbrdf.h"

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


