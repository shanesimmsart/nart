#pragma once

#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "sampling.h"

inline glm::vec3 Reflect(const glm::vec3& w1, const glm::vec3& w2)
{
    return (2.f * glm::dot(w1, w2) * w2) - w1;
}

// Dielectric Fresnel function
// Note: cosTheta is expected to be unsigned
float Fresnel(float eta_o, float eta_i, float cosTheta);

enum BSDFFlags
{
    SPECULAR = 1,
    GLOSSY = 2,
    DIFFUSE = 4,
    TRANSMISSIVE = 8
};

class BxDF
{
public:
    // BxDF - wo and wi must always be in local coord sys
    virtual glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime) = 0;

    // Sample BRDF, setting wi and pdf
    virtual glm::vec3 Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, float* alpha_i, bool use_alpha_prime) = 0;

    virtual float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime) = 0;

    virtual ~BxDF() {}

    uint8_t flags;

protected:
    BxDF() {}

    // alpha used when calculating BRDF
    float alpha = 0.f;
    // unadjusted alpha
    float alpha_0 = 0.f;
    // alpha adjusted via roughening over paths
    float alpha_prime = 0.f;
};

using BxDFPtr = std::unique_ptr<BxDF>;

class BSDF
{
public:
    BSDF(const glm::vec3& n, uint8_t numBxDFs);

    inline glm::vec3 ToLocal(const glm::vec3& v)
    {
        return glm::normalize(glm::vec3(glm::dot(v, n_t), glm::dot(v, n_b), glm::dot(v, n)));
    }

    inline glm::vec3 ToWorld(const glm::vec3& v)
    {
        return glm::normalize(v.x * n_t + v.y * n_b + v.z * n);
    }

    inline void AddBxDF(BxDFPtr&& bxdf)
    {
        bxdfs.push_back(std::move(bxdf));
    }

    // Sum BxDFs
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime);

    // Sample and sum BRDFs
    glm::vec3 Sample_f(const glm::vec3& wo, glm::vec3* wi, float sample1D, glm::vec2 sample, float* pdf, uint8_t* flags, bool use_alpha_prime, float* alpha_i = nullptr);

    // Average pdfs
    float Pdf(const glm::vec3& wo, const glm::vec3& wi, bool use_alpha_prime) const;

protected:
    // Coord sys in worldspace
    glm::vec3 n_t, n_b, n;
    uint8_t numBxDFs = 0;
    std::vector<BxDFPtr> bxdfs;
};


