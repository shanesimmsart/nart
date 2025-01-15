#include "../../include/nart/bxdfs/specularbrdf.h"

SpecularBRDF::SpecularBRDF(const glm::vec3& rho_s, float eta)
    : rho_s(rho_s), eta(eta) {
    flags = SPECULAR;
}

glm::vec3 SpecularBRDF::f(const glm::vec3& wo, const glm::vec3& wi,
                          bool use_alpha_prime) {
    // Probability of randomly sampling a delta function == 0
    return glm::vec3(0.f);
}

glm::vec3 SpecularBRDF::Sample_f(const glm::vec3& wo, glm::vec3& wi,
                                 float sample1D, glm::vec2 sample, float& pdf,
                                 uint8_t& flags, float* alpha_i,
                                 bool use_alpha_prime) {
    if (alpha_i != nullptr) *alpha_i = 0.f;
    flags = SPECULAR;

    wi = glm::vec3(-wo.x, -wo.y, wo.z);
    // Delta distribution, so PDF == 1 at this one sample point
    pdf = 1.f;

    // Mirror-reflection at grazing angles
    if (wi.z == 0.f) return glm::vec3(1.f);

    return (rho_s * Fresnel(1.f, eta, wi.z)) / glm::abs(wi.z);
}

float SpecularBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi,
                        bool use_alpha_prime) {
    return 0.f;
}
