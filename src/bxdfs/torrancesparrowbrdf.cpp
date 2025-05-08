#include "../../include/nart/bxdfs/torrancesparrowbrdf.h"

TorranceSparrowBRDF::TorranceSparrowBRDF(const glm::vec3& rho_s, float eta,
                                         float _alpha_0, float _alpha_prime)
    : rho_s(rho_s), eta(eta) {
    alpha = _alpha_0;
    alpha_0 = _alpha_0;
    alpha_prime = _alpha_prime;
    flags = GLOSSY;
}

float TorranceSparrowBRDF::Lambda(const glm::vec3& w) const {
    float sinTheta = glm::sqrt(1.f - (w.z * w.z));
    float tanTheta = (sinTheta / w.z);
    return (-1.f + glm::sqrt(1.f + (alpha * alpha * tanTheta * tanTheta))) *
           0.5f;
}

float TorranceSparrowBRDF::D(const glm::vec3 wh) const {
    float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
    float tanTheta = (sinTheta / wh.z);
    float tan2Theta = tanTheta * tanTheta;
    float theta = glm::asin(sinTheta);
    // I think I read somewhere that (x * x) * (x * x) is faster than
    // x * x * x * x
    return 1.f / ((glm::pi<float>() * alpha * alpha *
                   ((wh.z * wh.z) * (wh.z * wh.z))) *
                  (1.f + (tan2Theta / (alpha * alpha))) *
                  (1.f + (tan2Theta / (alpha * alpha))));
}

glm::vec3 TorranceSparrowBRDF::f(const glm::vec3& wo, const glm::vec3& wi,
                                 bool use_alpha_prime, float eta_outer) {
    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    if (wo.z < 0.f || wi.z < 0.f) return glm::vec3(0.f);

    glm::vec3 wh = wo + wi;
    wh = glm::normalize(wh);

    float g = G(wo, wi);
    float d = D(wh);
    float fr = Fresnel(eta_outer, eta, glm::dot(wh, wi));

    if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

    return (rho_s * g * d * fr) / (4.f * wo.z * wi.z);
}

glm::vec3 TorranceSparrowBRDF::Sample_f(const glm::vec3& wo, glm::vec3& wi,
                                        float sample1D, glm::vec2 sample,
                                        float& pdf, uint8_t& flags,
                                        float* alpha_i, bool use_alpha_prime, float eta_outer) {
    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    if (alpha_i != nullptr) *alpha_i = alpha;
    flags = SPECULAR;
    if (alpha > 0.001f) flags = GLOSSY;
    if (alpha >= 1.0f) flags = DIFFUSE;

    // Transform wo from ellipsoid to hemisphere
    glm::vec3 wo_h = glm::vec3(wo.x * alpha, wo.y * alpha, wo.z);
    wo_h = glm::normalize(wo_h);

    // Sample projection of hemisphere in wo
    // Build coord sys
    glm::vec3 T1 = glm::vec3(wo_h.y, -wo_h.x, 0.f);  // cross(wo, z)
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

    wi = Reflect(wo, wh);

    wi = glm::normalize(wi);

    pdf = Pdf(wo, wi, use_alpha_prime, eta_outer);

    return f(wo, wi, use_alpha_prime, eta_outer);
}

float TorranceSparrowBRDF::Get_eta() const { return eta; }

float TorranceSparrowBRDF::Pdf(const glm::vec3& wo, const glm::vec3& wi,
                               bool use_alpha_prime, float eta_outer) {
    if (use_alpha_prime)
        alpha = alpha_prime;
    else
        alpha = alpha_0;

    glm::vec3 wh = wo + wi;
    wh = glm::normalize(wh);

    if (wh.z < 0.f) return 0.f;

    float cosThetaH = glm::min(glm::dot(wo, wh), 1.f);
    float pdf = (D(wh) * glm::min(glm::dot(wo, wh), 1.f) * G1(wo)) / wo.z;
    return glm::max(0.f, pdf / (4.f * cosThetaH));
}
