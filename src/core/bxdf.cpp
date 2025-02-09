#include "../../include/nart/core/bxdf.h"

float Fresnel(float eta_o, float eta_i, float cosTheta) {
    float cosTheta_o = glm::min(glm::abs(cosTheta), 1.f);

    float sinTheta_o = glm::sqrt(1.f - (cosTheta_o * cosTheta_o));
    float sinTheta_i = (eta_o / eta_i) * sinTheta_o;

    // TIR
    if (sinTheta_i > 1.f) return 1.f;
    float cosTheta_i = glm::sqrt(1.f - (sinTheta_i * sinTheta_i));

    if (glm::abs(cosTheta_o + cosTheta_i) < 0.00001f) return 0.f;

    float fPara = ((eta_i * cosTheta_o) - (eta_o * cosTheta_i)) /
                  ((eta_i * cosTheta_o) + (eta_o * cosTheta_i));
    float fPerp = ((eta_o * cosTheta_o) - (eta_i * cosTheta_i)) /
                  ((eta_o * cosTheta_o) + (eta_i * cosTheta_i));
    return ((fPara * fPara) + (fPerp * fPerp)) * 0.5f;
}

BSDF::BSDF(const Intersection& isect, uint8_t numBxDFs)
    : n(isect.sn), numBxDFs(numBxDFs) {}

void BSDF::BuildCoordSys(const Intersection& isect, glm::vec3* nn) {
    // Build local coord sys from normal
    // Project dpds onto plane defined by n to get n_t
    n_t = glm::normalize(isect.dpds - glm::dot(isect.dpds, n) * n);

    // Get n_b
    n_b = glm::normalize(glm::cross(isect.sn, n_t));

    if (nn) {
        // Transform nn to world space
        n = glm::normalize(ToWorld(*nn));

        // Project dpds onto plane defined by nn to get new n_t
        n_t = glm::normalize(isect.dpds - glm::dot(isect.dpds, n) * n);

        // Get new n_b
        n_b = glm::normalize(glm::cross(isect.sn, n_t));
    }
}

glm::vec3 BSDF::f(const glm::vec3& wo, const glm::vec3& wi,
                  bool use_alpha_prime) {
    glm::vec3 f = glm::vec3(0.f);
    for (uint8_t i = 0; i < numBxDFs; ++i) {
        f += bxdfs[i]->f(wo, wi, use_alpha_prime);
    }
    return f;
}

glm::vec3 BSDF::Sample_f(const glm::vec3& wo, glm::vec3& wi, float sample1D,
                         glm::vec2 sample, float& pdf, uint8_t& flags,
                         bool use_alpha_prime, float* alpha_i) {
    // Choose a BxDF
    uint8_t bxdfIndex =
        static_cast<uint8_t>(sample1D * static_cast<float>(numBxDFs));

    // Remap sample to remove bias so we can reuse it
    sample1D = glm::fract(sample1D * static_cast<float>(numBxDFs));

    glm::vec3 f = bxdfs[bxdfIndex]->Sample_f(wo, wi, sample1D, sample, pdf,
                                             flags, alpha_i, use_alpha_prime);

    float numEvaluated = 1.f;

    if (!(flags & SPECULAR)) {
        for (uint8_t i = 0; i < numBxDFs; ++i) {
            if (i != bxdfIndex && !(bxdfs[i]->flags & SPECULAR)) {
                float bxdfPdf = bxdfs[i]->Pdf(wo, wi, use_alpha_prime);
                if (bxdfPdf > 0.f) {
                    pdf += bxdfs[i]->Pdf(wo, wi, use_alpha_prime);
                    f += bxdfs[i]->f(wo, wi, use_alpha_prime);
                    numEvaluated += 1.f;
                }
            }
        }
        pdf /= static_cast<float>(numBxDFs);  // numEvaluated;
    }

    return f;
}

float BSDF::Pdf(const glm::vec3& wo, const glm::vec3& wi,
                bool use_alpha_prime) const {
    float pdf = 0.f;

    for (uint8_t i = 0; i < numBxDFs; ++i) {
        pdf += bxdfs[i]->Pdf(wo, wi, use_alpha_prime);
    }

    return pdf / static_cast<float>(numBxDFs);
}

void BSDF::SetN(glm::vec3 *n) { nn = n; }

glm::vec3 BSDF::GetN(){
    return n;
}
