#include "reflection.h"

float Fresnel(float etaI, float etaT, float cosTheta)
{
    float cosThetaI = glm::min(cosTheta, 1.f);
    if (cosThetaI < 0.f) std::swap(etaI, etaT);

    float sinThetaI = glm::sqrt(1.f - (cosThetaI * cosThetaI));
    float sinThetaT = (etaI / etaT) * sinThetaI;
    // TIR
    if (sinThetaT > 1.f) return 1.f;
    float cosThetaT = glm::sqrt(1.f - (sinThetaT * sinThetaT));

    float fPara = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
    float fPerp = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
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