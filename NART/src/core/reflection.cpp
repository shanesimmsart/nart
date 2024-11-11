#include "reflection.h"

float Fresnel(float eta_o, float eta_i, float cosTheta)
{
    float cosTheta_o = glm::min(cosTheta, 1.f);
    if (cosTheta_o < 0.f)
    {
        cosTheta_o = -cosTheta_o;
        std::swap(eta_o, eta_i);
    }

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
        nt = glm::normalize(glm::vec3(n.z, 0.f, -n.x));
    }

    else
    {
        nt = glm::normalize(glm::vec3(0.f, n.z, -n.y));
    }
    nb = glm::normalize(glm::cross(n, nt));

    bxdfs.reserve(numBxDFs);
}