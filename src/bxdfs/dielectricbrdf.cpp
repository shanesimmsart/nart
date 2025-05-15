#include "../../include/nart/bxdfs/dielectricbrdf.h"

DielectricBRDF::DielectricBRDF(const glm::vec3& rho_s, const glm::vec3& tau,
                               float eta, float _alpha_0, float _alpha_prime)
    : rho_s(rho_s), tau(tau), eta(eta) {
    alpha = _alpha_0;
    alpha_0 = _alpha_0;
    alpha_prime = _alpha_prime;
    flags = GLOSSY;
}

float DielectricBRDF::Lambda(const glm::vec3& w) const {
    float sinTheta = glm::sqrt(1.f - (w.z * w.z));
    float tanTheta = (sinTheta / w.z);
    return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) *
           0.5f;
}

float DielectricBRDF::D(const glm::vec3 wh) const {
    if (wh.z == 0.f) return 0.f;
    float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
    float tanTheta = (sinTheta / wh.z);
    float tan2Theta = tanTheta * tanTheta;
    float theta = glm::asin(sinTheta);
    return 1.f / ((glm::pi<float>() * alpha * alpha *
                   ((wh.z * wh.z) * (wh.z * wh.z))) *
                  (1.f + (tan2Theta / (alpha * alpha))) *
                  (1.f + (tan2Theta / (alpha * alpha))));
}

glm::vec3 DielectricBRDF::f(const glm::vec3& wo, const glm::vec3& wi,
                            bool use_alpha_prime, float eta_outer) {
    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    float eta_o = eta_outer;
    float eta_i = eta;

    if (wo.z < 0.f) std::swap(eta_o, eta_i);

    if (wo.z * wi.z >= 0.f) {
        // Reflect
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);
        if (wh.z < 0.f) wh *= -1.f;

        float g = G(wo, wi);
        float d = D(wh);
        float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));

        if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

        return ((rho_s * g * d * Fr) / (4.f * wo.z * wi.z));
    }

    else {
        // Refract
        glm::vec3 wh = glm::normalize((eta_o * wo) + (eta_i * wi));
        if (wh.z < 0.f) wh *= -1.f;

        float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));
        if (Fr >= 1.f) {
            return glm::vec3(0.f);
        }

        float g = G(wo, wi);
        float d = D(wh);
        float wiDotWh = glm::dot(wi, wh);
        float woDotWh = glm::dot(wo, wh);

        float num = g * d * (1.f - Fr) * glm::abs(wiDotWh) *
                    glm::abs(woDotWh) * eta_o * eta_o;
        float denom = ((eta_i * wiDotWh) + (eta_o * woDotWh)) *
                      ((eta_i * wiDotWh) + (eta_o * woDotWh)) *
                      glm::abs(wo.z * wi.z);
        return (glm::vec3(num / denom)) * tau;
    }
}

glm::vec3 DielectricBRDF::Sample_f(const glm::vec3& wo, glm::vec3& wi,
                                   float sample1D, glm::vec2 sample, float& pdf,
                                   uint8_t& flags, float* alpha_i,
                                   bool use_alpha_prime, float eta_outer) {
    float eta_o = eta_outer;
    float eta_i = eta;

    if (eta_o == eta_i) {
        wi = -wo;
        pdf = 0.f;
        flags |= TRANSMISSIVE;
        return tau;
    }

    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    if (alpha_i != nullptr) *alpha_i = alpha;
    flags = SPECULAR;
    if (alpha > 0.0001f) flags = GLOSSY;
    if (alpha >= 1.0f) flags = DIFFUSE;

    // Transform wo from ellipsoid to hemisphere
    glm::vec3 wo_h = glm::vec3(wo.x * alpha, wo.y * alpha, wo.z);
    wo_h = glm::normalize(wo_h);
    if (wo.z < 0.f) wo_h *= -1.f;

    // Sample projection of hemisphere in wo
    // Build coord sys
    glm::vec3 T1;
    if (wo.x == 0.f && wo.y == 0.f)
        T1 = glm::vec3(1.f, 0.f, 0.f);
    else
        T1 = glm::vec3(wo_h.y, -wo_h.x, 0.f);  // cross(wo, z)
    T1 = glm::normalize(T1);
    glm::vec3 T2 = glm::cross(T1, wo_h);
    T2 = glm::normalize(T2);

    // Sample disk
    glm::vec2 visibleHSample = UniformSampleDisk(sample);
    // Convert to visible hemisphere coords
    float s = (1.f + wo_h.z) * 0.5f;
    visibleHSample.y =
        (s * visibleHSample.y) +
        ((1.f - s) * glm::sqrt(1.f - (visibleHSample.x * visibleHSample.x)));
    // Project onto hemisphere
    glm::vec3 wh =
        glm::vec3(glm::sqrt(1.f - (visibleHSample.x * visibleHSample.x) -
                            (visibleHSample.y * visibleHSample.y)),
                  visibleHSample.x, visibleHSample.y);
    // Convert to local coords
    wh = wh.x * wo_h + wh.y * T1 + wh.z * T2;

    // Transform wh to ellipsoid (Inverse transpose)
    wh = glm::vec3(wh.x * alpha, wh.y * alpha, wh.z);
    wh = glm::normalize(wh);

    // Refract
    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    float Fr = Fresnel(eta_o, eta_i, glm::abs(glm::dot(wh, wo)));

    if (sample1D < Fr) {
        // Reflect
        wi = Reflect(wo, wh);
        wi = glm::normalize(wi);

        pdf = Pdf(wo, wi, use_alpha_prime, eta_outer) * Fr;

        return f(wo, wi, use_alpha_prime, eta_outer);
    }

    // Refract
    float cosTheta_o = glm::min(1.f, glm::max(-1.f, glm::dot(wo, wh)));
    float sinTheta_o = glm::sqrt(1.f - (cosTheta_o * cosTheta_o));
    float sinTheta_i = ((eta_o / eta_i) * sinTheta_o);

    // TIR
    if (sinTheta_i >= 1.f) {
        // Reflect
        wi = Reflect(wo, wh);
        wi = glm::normalize(wi);

        pdf = Pdf(wo, wi, use_alpha_prime, eta_outer) * (1.f - Fr);

        return f(wo, wi, use_alpha_prime, eta_outer);
    }

    flags |= TRANSMISSIVE;

    glm::vec3 b = wh * cosTheta_o;
    glm::vec3 a = wo - b;
    glm::vec3 c = -a * (eta_o / eta_i);
    glm::vec3 d = -wh * glm::sqrt(1.f - (sinTheta_i * sinTheta_i));
    if (glm::dot(wo, wh) < 0.f) d *= -1.f;
    wi = glm::normalize(c + d);

    pdf = Pdf(wo, wi, use_alpha_prime, eta_outer) * (1.f - Fr);

    return f(wo, wi, use_alpha_prime, eta_outer);
}

float DielectricBRDF::Get_eta() const { return eta; }

float DielectricBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi,
                          bool use_alpha_prime, float eta_outer) {
    float eta_o = eta_outer;
    float eta_i = eta;

    if (eta_o == eta_i) {
        return 0.f;
    }

    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    if (wo.z * wi.z >= 0.f) {
        // Reflect
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        if (wh.z < 0.f) wh *= -1.f;

        float cosThetaH = glm::abs(glm::min(glm::dot(wo, wh), 1.f));
        float pdf = (D(wh) * glm::min(glm::dot(wo, wh), 1.f) * G1(wo)) / wo.z;
        return glm::max(0.f, pdf / (4.f * cosThetaH));
    }

    // Refract
    if (wo.z < 0.f) std::swap(eta_o, eta_i);
    glm::vec3 wh = glm::normalize((eta_o * wo) + (eta_i * wi));
    if (wh.z < 0.f) wh *= -1.f;

    float pdf = (D(wh) * glm::min(glm::abs(glm::dot(wo, wh)), 1.f) * G1(wo)) /
                glm::abs(wo.z);
    float dotWiWh = glm::dot(wi, wh);
    float dotWoWh = glm::dot(wo, wh);
    float denom = (eta_i * dotWiWh + eta_o * dotWoWh);
    float JDet = (std::abs(dotWiWh) * eta_i * eta_i) / (denom * denom);
    return pdf * JDet;
}
