#include "../../include/nart/bxdfs/speculardielectricbrdf.h"

SpecularDielectricBRDF::SpecularDielectricBRDF(const glm::vec3& rho_s,
                                               const glm::vec3& tau, float eta)
    : rho_s(rho_s), tau(tau), eta(eta) {
    flags = SPECULAR;
}

glm::vec3 SpecularDielectricBRDF::f(const glm::vec3& wo, const glm::vec3& wi,
                                    bool use_alpha_prime, float eta_outer) {
    // Probability of randomly sampling a delta function == 0
    return glm::vec3(0.f);
}

glm::vec3 SpecularDielectricBRDF::Sample_f(const glm::vec3& wo, glm::vec3& wi,
                                           float sample1D, glm::vec2 sample,
                                           float& pdf, uint8_t& flags,
                                           float* alpha_i,
                                           bool use_alpha_prime, float eta_outer) {
    float eta_o = eta_outer;
    float eta_i = eta;

    if (eta_o == eta_i) {
        wi = -wo;
        pdf = 0.f;
        flags |= TRANSMISSIVE;
        return tau;
    }
    
    if (alpha_i != nullptr) *alpha_i = 0.f;
    flags = SPECULAR;

    // Check if inside or outside of medium
    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    float Fr = Fresnel(eta_o, eta_i, glm::abs(wo.z));

    if (sample.x < Fr) {
        // Delta distribution, so PDF == 1 at this one sample point
        pdf = Fr;

        // Reflect
        wi = glm::vec3(-wo.x, -wo.y, wo.z);

        // Mirror-reflection at grazing angles
        if (wi.z == 0.f) return glm::vec3(1.f);

        // Check hemisphere with gn
        return glm::vec3(Fr / glm::abs(wi.z)) * rho_s;
    }

    else {
        // Delta distribution, so PDF == 1 at this one sample point
        pdf = 1.f - Fr;

        // Refract
        float sinTheta_o = glm::sqrt(1.f - (wo.z * wo.z));
        float sinTheta_i = ((eta_o / eta_i) * sinTheta_o);

        // TIR
        if (sinTheta_i >= 1.f) {
            wi = glm::vec3(-wo.x, -wo.y, wo.z);
            return glm::vec3(1.f) * rho_s;
        }

        flags |= TRANSMISSIVE;

        glm::vec3 n(0.f, 0.f, 1.f);
        glm::vec3 b = n * wo.z;
        glm::vec3 a = wo - b;
        glm::vec3 c = -a * (eta_o / eta_i);
        glm::vec3 d = -n * glm::sqrt(1.f - (sinTheta_i * sinTheta_i));
        if (wo.z < 0.f) d *= -1.f;
        wi = glm::normalize(c + d);

        // Check hemisphere with gn
        glm::vec3 f = glm::vec3(
            ((1.f - Fr)) / glm::abs(wi.z)); // (eta_o / eta_i) * (eta_o / eta_i) * 
        f *= tau;

        return f;
    }
}

float SpecularDielectricBRDF::Get_eta() const { return eta; }

float SpecularDielectricBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi,
                                  bool use_alpha_prime, float eta_outer) {
    return 0.f;
}
