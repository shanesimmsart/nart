#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <queue>
#include <tbb/task_group.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>
#include <glm/common.hpp>
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/trigonometric.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/constants.hpp>
#include <nlohmann/json.hpp>

#define Infinity std::numeric_limits<float>::infinity()
#define FilterTableResolution 64

#define DEBUG_BUCKET 0
#define DEBUG_BUCKET_X 53
#define DEBUG_BUCKET_Y 21

#define BSDF_SAMPLING 1
#define LIGHT_SAMPLING 1


glm::vec2 UniformSampleDisk(glm::vec2 sample)
{
    float r = glm::sqrt(sample.x);
    float theta = sample.y * glm::two_pi<float>();

    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);

    float x = r * cosTheta;
    float y = r * sinTheta;

    return glm::vec2(x, y);
}

glm::vec2 UniformSampleRing(glm::vec2 sample, float* pdf, float innerRadius)
{
    // float r = glm::sqrt(sample.x * (1.f - r2) * (1.f + r2));
    // float r = glm::sqrt(sample.x * (1.f - (r2 * r2)));
    float r = glm::sqrt(glm::mix(innerRadius, 1.f, sample.x));
    float theta = sample.y * glm::two_pi<float>();

    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);

    float x = r * cosTheta;
    float y = r * sinTheta;

    *pdf = 1.f / (glm::pi<float>() * (1.f - innerRadius));

    return glm::vec2(x, y);
}

glm::vec3 CosineSampleHemisphere(glm::vec2 sample, float* pdf)
{
    // Using Malley's method, we can uniformly sample a disk, then project
    // the sample points upwards onto the hemisphere to achieve a
    // cosine-weighted distribution
    glm::vec2 diskSample = UniformSampleDisk(sample);
    float z = glm::sqrt(1.f - (diskSample.x * diskSample.x + diskSample.y * diskSample.y));

    *pdf = z * glm::one_over_pi<float>();

    return glm::vec3(diskSample, z);
}

// Information about the render settings to be used
struct RenderInfo
{
    RenderInfo(uint32_t imageWidth, uint32_t imageHeight, uint32_t bucketSize, uint32_t spp, float filterWidth) :
        imageWidth(imageWidth), imageHeight(imageHeight), bucketSize(bucketSize), spp(spp), filterWidth(filterWidth)
    {
        filterBounds = static_cast<uint32_t>(glm::ceil(filterWidth));
        tileSize = bucketSize + (filterBounds * 2);
        totalWidth = imageWidth + (filterBounds * 2);
        totalHeight = imageHeight + (filterBounds * 2);
    }

    const uint32_t imageWidth;
    const uint32_t imageHeight;
    const uint32_t bucketSize;
    const uint32_t spp;
    const float filterWidth;

    // Discrete bounds of filter
    uint32_t filterBounds;
    // Size of each tile, including border added to include filter width
    uint32_t tileSize;
    // Total image dimensions, including border added to include filter width
    uint32_t totalWidth;
    uint32_t totalHeight;
};

class Ray
{
public:
    Ray(glm::vec3 o, glm::vec3 d) : o(o), d(d) {}

    glm::vec3 o = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 d = glm::vec3(0.f, 1.f, 0.f);
};

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

enum BSDFFlags
{
    SPECULAR = 1,
    GLOSSY = 2,
    DIFFUSE = 4
};

class BxDF
{
public:
    // BxDF - wo and wi must always be in local coord sys
    virtual glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi) = 0;

    // Sample BRDF, setting wi and pdf
    virtual glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags) = 0;

    virtual float Pdf(const glm::vec3& wo, const glm::vec3& wi) = 0;

    BSDFFlags flags;

protected:
    BxDF() {}
};

class BSDF
{
public:
    BSDF(glm::vec3 n, uint8_t numBxDFs) : n(n), numBxDFs(numBxDFs)
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

    glm::vec3 ToLocal(glm::vec3 v)
    {
        return glm::normalize(glm::vec3(glm::dot(v, nt), glm::dot(v, nb), glm::dot(v, n)));
    }

    glm::vec3 ToWorld(glm::vec3 v)
    {
        return glm::normalize(v.x * nt + v.y * nb + v.z * n);
    }

    void AddBxDF(std::shared_ptr<BxDF> bxdf)
    {
        bxdfs.emplace_back(bxdf);
    }

    // Sum BxDFs
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        glm::vec3 f = glm::vec3(0.f);
        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            f += bxdfs[i]->f(wo, wi);
        }
        return f;
    }

    // Sample and sum BRDFs
    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        // Choose a BxDF
        uint8_t bxdfIndex = static_cast<uint8_t>(sample.x * static_cast<float>(numBxDFs));

        // Remap sample to remove bias
        sample.x = glm::fract(sample.x * static_cast<float>(numBxDFs));

        glm::vec3 f = bxdfs[bxdfIndex]->Sample_f(wo, wi, sample, pdf, flags);

        float numEvaluated = 1.f;

        if (!(*flags & SPECULAR))
        {
            for (uint8_t i = 0; i < numBxDFs; ++i)
            {
                if (i != bxdfIndex && !(bxdfs[i]->flags & SPECULAR))
                {
                    float bxdfPdf = bxdfs[i]->Pdf(wo, *wi);
                    if (bxdfPdf > 0.f)
                    {
                        *pdf += bxdfs[i]->Pdf(wo, *wi);
                        f += bxdfs[i]->f(wo, *wi);
                        numEvaluated += 1.f;
                    }
                }
            }
            *pdf /= static_cast<float>(numBxDFs); //numEvaluated;
        }

        return f;
    }

    // Average pdfs
    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        float pdf = 0.f;

        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            pdf += bxdfs[i]->Pdf(wo, wi);
        }

        return pdf / static_cast<float>(numBxDFs);
    }

protected:
    // Coord sys in worldspace
    glm::vec3 nt, nb, n;
    uint8_t numBxDFs = 0;
    std::vector<std::shared_ptr<BxDF>> bxdfs;
};

class LambertBRDF : public BxDF
{
public:
    LambertBRDF(glm::vec3 rho) : rho(rho)
    {
        flags = DIFFUSE;
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        return rho * glm::one_over_pi<float>();
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = DIFFUSE;
        *wi = CosineSampleHemisphere(sample, pdf);
        return f(wo, *wi);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        return wi.z * glm::one_over_pi<float>();
    }

private:
    // Albedo
    glm::vec3 rho;
};

class SpecularBRDF : public BxDF
{
public:
    SpecularBRDF(glm::vec3 R, float eta) : R(R), eta(eta)
    {
        flags = SPECULAR;
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        // Probability of randomly sampling a delta function == 0
        return glm::vec3(0.f);
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = SPECULAR;

        *wi = glm::vec3(-wo.x, -wo.y, wo.z);
        // Delta distribution, so PDF == 1 at this one sample point
        *pdf = 1.f;

        // Mirror-reflection at grazing angles
        if (wi->z == 0.f) return glm::vec3(1.f);

        return (R * Fresnel(1.f, eta, wi->z)) / glm::abs(wi->z);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        return 0.f;
    }

private:
    // Reflectance
    glm::vec3 R;
    // IOR
    float eta;
};

glm::vec3 Reflect(const glm::vec3& w1, const glm::vec3& w2)
{
    return (2.f * glm::dot(w1, w2) * w2) - w1;
}

float RoughnessToAlpha(float roughness)
{
    return 1.62142f * glm::sqrt(roughness);
}

float AlphaToRoughness(float alpha)
{
    return (alpha * alpha) * 0.38037235782f;
}

class TorranceSparrowBRDF : public BxDF
{
public:
    TorranceSparrowBRDF(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
    {
        flags = GLOSSY;
    }

    // Masking-shadowing Function (Smith) (GGX)
    float Lambda(const glm::vec3& w)    
    {
        float sinTheta = glm::sqrt(1.f - (w.z * w.z));
        float tanTheta = (sinTheta / w.z);
        float a = 1.f / (alpha * tanTheta);
        return (std::erf(a) - 1.f + (glm::exp(-a * a) / (a * glm::sqrt(glm::pi<float>())))) * 0.5f;
    }

    float G(const glm::vec3& wo, const glm::vec3& wi)
    {
        return 1.f / (1.f + Lambda(wo) + Lambda(wi));
    }

    float G1(const glm::vec3& w)
    {
        return 1.f / (1.f + Lambda(w));
    }

    // Normal Distribution Function (GGX)
    float D(const glm::vec3 wh)
    {
        float sinTheta = glm::sqrt(1.f - (wh.z * wh.z));
        float tanTheta = (sinTheta / wh.z);
        float tan2Theta = tanTheta * tanTheta;
        // I think I read somewhere that (x * x) * (x * x) is faster than
        // x * x * x * x
        float theta = glm::asin(sinTheta);
        return (glm::exp(-tan2Theta / (alpha * alpha))) / (glm::pi<float>() * alpha * alpha * ((wh.z * wh.z) * (wh.z * wh.z)));
    }

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        if (wo.z < 0.f || wi.z < 0.f) return glm::vec3(0.f);

        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        float g = G(wo, wi);
        float d = D(wh);
        float fr = Fresnel(1.f, eta, glm::dot(wh, wi));

        if (wo.z * wi.z == 0.f) return glm::vec3(0.f);

        return (R * g * d * fr) / (4.f * wo.z * wi.z);
        // return (R  * D(wh)) / (4.f * wo.z * wi.z);
    }

    glm::vec3 Sample_f(glm::vec3 wo, glm::vec3* wi, glm::vec2 sample, float* pdf, BSDFFlags* flags)
    {
        *flags = SPECULAR;
        if (alpha > 0.001f) *flags = GLOSSY;
        // TODO: PBRT uses this as roughness == 1
        // Where does this come from??
        if (alpha >= 1.62142f) *flags = DIFFUSE;

        // Sample distribution to get wh
        float tanTheta = glm::sqrt(-alpha * alpha * glm::log(1.f - sample.x));
        float phi = 2.f * glm::pi<float>() * sample.y;

        // TODO: There is definitely a better way to do this...
        float theta = glm::atan(tanTheta);
        float x = glm::cos(phi) * glm::sin(theta);
        float y = glm::sin(phi) * glm::sin(theta);
        float z = glm::cos(theta);

        if (z < 0.f)
        {
            *pdf = 0.f;
            return glm::vec3(0.f);
        }

        glm::vec3 wh(x, y, z);

        *wi = Reflect(wo, wh);

        *wi = glm::normalize(*wi);

        *pdf = Pdf(wo, *wi);

        return f(wo, *wi);
    }

    float Pdf(const glm::vec3& wo, const glm::vec3& wi)
    {
        glm::vec3 wh = wo + wi;
        wh = glm::normalize(wh);

        if (wh.z < 0.f) return 0.f;

        // Need to convert 1/dwh to 1/dwi
        float cosThetaH = glm::min(glm::dot(wh, wi), 1.f);
        float d = D(wh);
        return (wh.z * d) / (4.f * cosThetaH);
    }

private:
    // Reflectance
    glm::vec3 R;
    // IOR
    float eta;
    // Roughness
    float alpha;
};

class Material
{
public:
    virtual BSDF CreateBSDF(glm::vec3 n, float roughnessOffset) = 0;

protected:
    Material() {};
};

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial(glm::vec3 rho) : rho(rho) // your boat...
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);
        bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));
        return bsdf;
    }

private:
    // Eventually, material inputs will be replaced with patterns, that vary depending on UVs, position, etc.
    glm::vec3 rho;
};

class SpecularMaterial : public Material
{
public:
    SpecularMaterial(glm::vec3 R, float eta) : R(R), eta(eta)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);

        float alpha = glm::min(roughnessOffset, 1.f);

        if (alpha > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alpha));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    glm::vec3 R;
    float eta;
};

class GlossyDielectricMaterial : public Material
{
public:
    GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha) : R(R), eta(eta), alpha(alpha)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 1);

        float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

        if (alphaAdjusted > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    float alpha;
    glm::vec3 R;
    float eta;
};

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha) : rho(rho), R(R), eta(eta), alpha(alpha)
    {}

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset)
    {
        BSDF bsdf(n, 2);

        bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));

        float alphaAdjusted = glm::min(alpha + roughnessOffset, 1.f);

        if (alphaAdjusted > 0.001f) {
            bsdf.AddBxDF(std::make_shared<TorranceSparrowBRDF>(R, eta, alphaAdjusted));
        }

        else bsdf.AddBxDF(std::make_shared<SpecularBRDF>(R, eta));

        return bsdf;
    }

private:
    glm::vec3 rho;
    float alpha;
    glm::vec3 R;
    float eta;
};

struct Intersection
{
    std::shared_ptr<Material> material;
    // Barycentric coords
    float u, v;
    // Point, geometric normal, shading normal
    glm::vec3 p, gn, sn;

    // Upper and lower bounds of ray to be considered
    float tMin = 0.f;
    float tMax = Infinity;

    // Flag for area light intersection
    bool isLight = false;
};

class Light
{
public:
    Light(glm::vec3 Le, float intensity, glm::mat4 LightToWorld) : Le(Le), intensity(intensity), LightToWorld(LightToWorld) {}

    // Return incident radiance from light, set wi
    virtual glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) = 0;

    // Sample light, return incident radiance, set wi, return pdf
    virtual glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf) = 0;

    virtual float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) = 0;

    // Does light have a delta distribution?
    bool isDelta = false;

protected:
    // Radiance emitted
    glm::vec3 Le;
    float intensity;
    glm::mat4 LightToWorld;
};

class DistantLight : public Light
{
public:
    DistantLight(glm::vec3 Le, float intensity, glm::mat4 LightToWorld) : Light(Le, intensity, LightToWorld)
    {
        isDelta = true;
        // Pointing down by default
        direction = glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld;
    }

    glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        return glm::vec3(0.f);
    }

    glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf)
    {
        *wi = -direction;
        *pdf = 1.f;
        return Le * intensity;
    }

    float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        return 0.f;
    }

    glm::vec3 direction;
};

class DiskLight : public Light
{
public:
    DiskLight(float radius, glm::vec3 Le, float intensity, glm::mat4 LightToWorld) : Light(Le, intensity, LightToWorld), radius(radius)
    {
        isDelta = false;
    }

    glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        if (Pdf(lightIsect, p, wi) > 0.f) return Le * intensity;

        else return glm::vec3(0.f);
    }

    glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf)
    {
        // Sample disk
        glm::vec4 diskSample = glm::vec4(UniformSampleDisk(sample) * radius, 0.f, 1.f);

        // Transform disk sample to world space
        diskSample = diskSample * LightToWorld;
        glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

        // Set w
        *wi = glm::vec3(diskSample) - p;
        float distPToSample = glm::sqrt(wi->x * wi->x + wi->y * wi->y + wi->z * wi->z);
        *wi = glm::normalize(*wi);

        // Calculate pdf with respect to disk area
        *pdf = 1.f / (glm::pi<float>() * radius * radius);

        float wiDotN = glm::dot(-*wi, n);
        if (wiDotN <= 0.f)
        {
            *pdf = 0.f;
            return glm::vec3(0.f);
        }
        // Calculate pdf projected onto hemisphere around p
        *pdf = *pdf * ((distPToSample * distPToSample) / wiDotN);

        lightIsect->p = diskSample;
        lightIsect->tMax = distPToSample;

        return Le * intensity;
    }

    float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        glm::vec3 diskCenter = glm::vec3(glm::vec4(0.f, 0.f, 0.f, 1.f) * LightToWorld);
        glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

        if (glm::dot(wi, n) >= 0.f) return 0.f;

        // Distance of plane from origin
        float D = glm::dot(diskCenter, n);

        // Intersect plane
        float t = (D - glm::dot(p, n)) / glm::dot(wi, n);
        if (t < 0.f) return 0.f;

        glm::vec3 pHit = p + t * wi;

        glm::vec3 diskCenterToPHit = pHit - diskCenter;

        // If not hit, pdf == 0
        float dist = diskCenterToPHit.x * diskCenterToPHit.x + diskCenterToPHit.y * diskCenterToPHit.y + diskCenterToPHit.z * diskCenterToPHit.z;
        if (dist > radius * radius) return 0.f;

        // If hit, calculate pdf projected onto hemisphere around p
        float pdf = 1.f / (glm::pi<float>() * radius * radius);
        pdf = pdf * ((t * t) / glm::dot(-wi, n));

        // TODO: Figure out why this is negative??
        lightIsect->tMax = t;

        return pdf;
    }

    float radius = 1.f;
};

class RingLight : public Light
{
public:
    RingLight(float radius, float innerRadius, glm::vec3 Le, float intensity, glm::mat4 LightToWorld) : Light(Le, intensity, LightToWorld), radius(radius), innerRadius(innerRadius)
    {
        isDelta = false;
    }

    glm::vec3 Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        if (Pdf(lightIsect, p, wi) > 0.f) return Le * intensity;

        else return glm::vec3(0.f);
    }

    glm::vec3 Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf)
    {
        // Sample disk
        glm::vec4 ringSample = glm::vec4(UniformSampleRing(sample, pdf, innerRadius / radius) * radius, 0.f, 1.f);

        // Transform disk sample to world space
        ringSample = ringSample * LightToWorld;
        glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

        // Set wi
        *wi = glm::vec3(ringSample) - p;
        *wi = glm::vec3(ringSample) - p;
        float distPToSample = glm::sqrt(wi->x * wi->x + wi->y * wi->y + wi->z * wi->z);
        *wi = glm::normalize(*wi);

        // Calculate pdf with respect to disk area
        *pdf /= (glm::pi<float>() * radius * radius);

        // Calculate pdf projected onto hemisphere around p
        float wiDotN = glm::dot(-*wi, n);
        if (wiDotN <= 0.f)
        {
            *pdf = 0.f;
            return glm::vec3(0.f);
        }
        // Calculate pdf projected onto hemisphere around p
        *pdf = *pdf * ((distPToSample * distPToSample) / wiDotN);

        lightIsect->p = ringSample;
        lightIsect->tMax = distPToSample;

        return Le * intensity;
    }

    float Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi)
    {
        glm::vec3 ringCenter = glm::vec3(glm::vec4(0.f, 0.f, 0.f, 1.f) * LightToWorld);
        glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

        if (glm::dot(wi, n) >= 0.f) return 0.f;

        // Distance of plane from origin
        float D = glm::dot(ringCenter, n);

        // Intersect plane
        float t = (D - glm::dot(p, n)) / glm::dot(wi, n);
        if (t < 0.f) return 0.f;

        glm::vec3 pHit = p + t * wi;

        glm::vec3 diskCenterToPHit = pHit - ringCenter;

        // If not hit, pdf == 0
        float dist = diskCenterToPHit.x * diskCenterToPHit.x + diskCenterToPHit.y * diskCenterToPHit.y + diskCenterToPHit.z * diskCenterToPHit.z;
        if (dist > radius * radius) return 0.f;
        if (dist < innerRadius * innerRadius) return 0.f;

        // If hit, calculate pdf projected onto hemisphere around p
        float pdf = 1.f / (glm::pi<float>() * (1.f - ((innerRadius * innerRadius) / (radius * radius))) * radius * radius);
        pdf = pdf * ((t * t) / glm::dot(-wi, n));

        lightIsect->tMax = t;

        return pdf;
    }

    float innerRadius = 0.f;
    float radius = 1.f;
};

class Triangle
{
public:
    Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2) :
        v0(v0), v1(v1), v2(v2), n0(n0), n1(n1), n2(n2) {}
    bool Intersect(const Ray& ray, Intersection* isect)
    {
        // Un-normalised direction of surface
        glm::vec3 n = glm::cross((v2 - v1), (v0 - v1));
        // Length of cross product == area of triangle bounded by vectors * 2
        float area = glm::length(n) / 2.f;

        // Back-face culling
        if (glm::dot(n, ray.d) >= 0.f) return false;

        // Find t of intersection with plane
        // Note: n doesn't need to be normalised yet as denominator divides by n's length
        float t = (glm::dot(v0, n) - glm::dot(ray.o, n)) / glm::dot(ray.d, n);

        if (t <= isect->tMin || t >= isect->tMax) return false;

        glm::vec3 p = ray.o + t * ray.d;

        // Check if p is on the correct side of each vector, and compute barycentric coords
        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0p = p - v0;
        if (glm::dot((glm::cross(v0v1, v0p)), n) < 0.f) return false;

        glm::vec3 v1v2 = v2 - v1;
        glm::vec3 v1p = p - v1;
        glm::vec3 c1 = glm::cross(v1v2, v1p);
        if (glm::dot(c1, n) < 0.f) return false;
        float u = glm::length(c1) / 2.f / area;

        glm::vec3 v2v0 = v0 - v2;
        glm::vec3 v2p = p - v2;
        glm::vec3 c2 = glm::cross(v2v0, v2p);
        if (glm::dot(c2, n) < 0.f) return false;

        isect->tMax = t;
        isect->u = u;
        isect->v = glm::length(c2) / 2.f / area;
        isect->gn = glm::normalize(n);
        isect->sn = n0 * isect->u + n1 * isect->v + n2 * (1 - isect->u - isect->v);
        isect->p = p;

        return true;
    }

    const glm::vec3 v0, v1, v2;
    const glm::vec3 n0, n1, n2;
};

class TriMesh
{
public:
    TriMesh(std::shared_ptr<Material> material) : material(material) {}

    // Triangles are passed directly, outside of the ctor, to reduce memory allocations
    std::vector<Triangle> triangles;
    std::shared_ptr<Material> material;
};

struct BoundingVolume
{
    BoundingVolume()
    {
        // Bounding volumes are defined by slabs along the cardinal direction vectors
        // and diagonally oriented vectors
        float OneOverSqrt3 = 1.f / glm::sqrt(3);
        normals[0] = glm::vec3(1, 0, 0);
        normals[1] = glm::vec3(0, 1, 0);
        normals[2] = glm::vec3(0, 0, 1);
        normals[3] = glm::vec3(OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[4] = glm::vec3(OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);
        normals[5] = glm::vec3(-OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[6] = glm::vec3(-OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);

        for (uint8_t i = 0; i < 7; ++i)
        {
            boundsMin[i] = Infinity;
            boundsMax[i] = -Infinity;
        }
    }

    // Extend bounding volume by another to ensure it is fully contained
    void ExtendBy(const BoundingVolume& bv)
    {
        for (uint8_t i = 0; i < 7; ++i)
        {
            boundsMin[i] = glm::min(boundsMin[i], bv.boundsMin[i]);
            boundsMax[i] = glm::max(boundsMax[i], bv.boundsMax[i]);
        }
    }

    bool Intersect(const Ray& ray, Intersection* isect)
    {
        float tMin = -Infinity;
        float tMax = Infinity;
        float dotLow = Infinity;

        for (uint8_t i = 0; i < 7; ++i)
        {
            // Equation of a plane:
            // (o+td)dotN - bound = 0
            // OdotN + tDdotN = bound
            // tDdotN = bound - OdotN
            // t = (bound - OdotN) / DdotN

            float slabMin = (boundsMin[i] - glm::dot(ray.o, normals[i])) / glm::dot(ray.d, normals[i]);
            float slabMax = (boundsMax[i] - glm::dot(ray.o, normals[i])) / glm::dot(ray.d, normals[i]);

            if (slabMin > slabMax)
            {
                std::swap(slabMin, slabMax);
            }

            if (slabMin > tMax || tMin > slabMax)
            {
                return false;
            }

            if (glm::dot(ray.o + (slabMin * ray.d), normals[i]) < dotLow)
            {
                dotLow = glm::dot(ray.o + (slabMin * ray.d), normals[i]);
            }

            tMin = glm::max(tMin, slabMin);
            tMax = glm::min(tMax, slabMax);
        }

        isect->tMin = tMin;
        isect->tMax = tMax;

        return true;
    }

    glm::vec3 normals[7];
    float boundsMin[7];
    float boundsMax[7];
};

// A chunk of triangles to be placed inside of a bounding volume
class Chunk
{
    // Stores a mesh pointer and an index into one of it's triangles
    struct TriangleIndex
    {
        TriangleIndex(std::shared_ptr<TriMesh> mesh, const uint32_t& index) : mesh(mesh), index(index) {}
        std::shared_ptr<TriMesh> mesh;
        const uint32_t index;
    };

public:
    void Insert(std::shared_ptr<TriMesh> mesh, const uint32_t& index)
    {
        triangleIndices.push_back(TriangleIndex(mesh, index));
    }

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        bool hit = false;
        for (uint32_t i = 0; i < triangleIndices.size(); ++i)
        {
            if (triangleIndices[i].mesh->triangles[triangleIndices[i].index].Intersect(ray, isect))
            {
                hit = true;
                isect->material = triangleIndices[i].mesh->material;
            }
        }
        return hit;
    }

    void CalculateBounds()
    {
        float OneOverSqrt3 = 1.f / glm::sqrt(3);

        glm::vec3 normals[7];
        normals[0] = glm::vec3(1, 0, 0);
        normals[1] = glm::vec3(0, 1, 0);
        normals[2] = glm::vec3(0, 0, 1);
        normals[3] = glm::vec3(OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[4] = glm::vec3(OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);
        normals[5] = glm::vec3(-OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[6] = glm::vec3(-OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);

        for (TriangleIndex triIndex : triangleIndices)
        {
            Triangle triangle = triIndex.mesh->triangles[triIndex.index];

            bboxMin = glm::min(bboxMin, triangle.v0);
            bboxMax = glm::max(bboxMax, triangle.v0);

            bboxMin = glm::min(bboxMin, triangle.v1);
            bboxMax = glm::max(bboxMax, triangle.v1);

            bboxMin = glm::min(bboxMin, triangle.v2);
            bboxMax = glm::max(bboxMax, triangle.v2);

            for (uint8_t i = 0; i < 7; ++i)
            {
                float v0DotNorm = glm::dot(triangle.v0, normals[i]);
                bv.boundsMin[i] = glm::min(v0DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v0DotNorm, bv.boundsMax[i]);

                float v1DotNorm = glm::dot(triangle.v1, normals[i]);
                bv.boundsMin[i] = glm::min(v1DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v1DotNorm, bv.boundsMax[i]);

                float v2DotNorm = glm::dot(triangle.v2, normals[i]);
                bv.boundsMin[i] = glm::min(v2DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v2DotNorm, bv.boundsMax[i]);
            }
        }
    }

    // Bounding box (used for octree position)
    glm::vec3 bboxMin = glm::vec3(Infinity);
    glm::vec3 bboxMax = glm::vec3(-Infinity);
    // Bounding volume (used for intersecting)
    BoundingVolume bv = BoundingVolume();

private:
    std::vector<TriangleIndex> triangleIndices;
};

// BVH uses an octree structure for containing bounding volumes
class Octree
{
    struct OctreeNode
    {
        OctreeNode(glm::vec3 nodeMin, glm::vec3 nodeMax) : nodeMin(nodeMin), nodeMax(nodeMax) {}

        std::vector<std::shared_ptr<Chunk>> chunks;
        std::shared_ptr<OctreeNode> children[8] = { nullptr };
        bool isLeaf = true;

        const glm::vec3 nodeMin;
        const glm::vec3 nodeMax;

        BoundingVolume bv = BoundingVolume();
    };

public:
    Octree(glm::vec3 sceneMin, glm::vec3 sceneMax) : sceneMin(sceneMin), sceneMax(sceneMax)
    {
        root = std::make_shared<OctreeNode>(sceneMin, sceneMax);
    }

    void Insert(std::shared_ptr<Chunk> chunk)
    {
        chunk->CalculateBounds();
        InsertChunkIntoNode(chunk, root);
    }

    void Build()
    {
        BuildBoundingVolumes(root);
    }

    bool Intersect(const Ray& ray, Intersection* isect)
    {
        // Using priority queue to keep track of closest node intersection
        std::priority_queue<std::pair<float, std::shared_ptr<OctreeNode>>, std::vector<std::pair<float, std::shared_ptr<OctreeNode>>>, std::greater<std::pair<float, std::shared_ptr<OctreeNode>>>> queue;
        queue.push(std::make_pair(Infinity, root));

        bool hit = false;

        while (!queue.empty())
        {
            std::shared_ptr<OctreeNode> topNode = queue.top().second;
            queue.pop();
            Intersection triIsect;
            for (uint8_t i = 0; i < 8; ++i)
            {
                std::shared_ptr<OctreeNode> node = topNode->children[i];
                if (node)
                {
                    Intersection bvIsect;
                    if (node->bv.Intersect(ray, &bvIsect))
                    {
                        queue.push(std::make_pair(bvIsect.tMin, topNode->children[i]));

                        if (node->isLeaf)
                        {
                            for (uint32_t j = 0; j < node->chunks.size(); ++j)
                            {
                                if (node->chunks[j]->Intersect(ray, &triIsect))
                                {
                                    if (triIsect.tMax < isect->tMax) hit = true;
                                }
                            }
                        }
                    }
                }
            }
            // If the nearest triangle is closer than the nearest bounding volume in the queue, we don't need to continue
            if (!queue.empty() && isect->tMax < queue.top().first) break;
            if (isect->tMax > triIsect.tMax) *isect = triIsect;
        }

        return hit;
    }

private:
    void InsertChunkIntoNode(std::shared_ptr<Chunk> chunk, std::shared_ptr<OctreeNode> node, uint8_t depth = 0)
    {
        if (node->isLeaf)
        {
            if (node->chunks.empty() || depth >= maxDepth)
            {
                node->chunks.push_back(chunk);
            }

            else
            {
                node->isLeaf = false;

                for (uint32_t i = 0; i < node->chunks.size(); ++i)
                {
                    InsertChunkIntoNode(node->chunks.back(), node, depth++);
                    node->chunks.pop_back();
                }

                InsertChunkIntoNode(chunk, node, depth++);
            }
        }

        else
        {
            // Figure out child node index based on chunk centroid and node centroid
            uint8_t nodeIndex = 0;

            glm::vec3 chunkCentroid = chunk->bboxMin + ((chunk->bboxMax - chunk->bboxMin) * 0.5f);
            glm::vec3 nodeCentroid = node->nodeMin + ((node->nodeMax - node->nodeMin) * 0.5f);

            if (chunkCentroid.x > nodeCentroid.x) nodeIndex = nodeIndex | 1;
            if (chunkCentroid.y > nodeCentroid.y) nodeIndex = nodeIndex | 2;
            if (chunkCentroid.z > nodeCentroid.z) nodeIndex = nodeIndex | 4;

            if (!node->children[nodeIndex])
            {
                glm::vec3 childNodeSize = (node->nodeMax - node->nodeMin) * 0.5f;
                glm::vec3 childNodeMin = nodeCentroid;
                glm::vec3 childNodeMax = nodeCentroid;

                // Compute child bounds
                if (nodeIndex & 1) childNodeMax.x += childNodeSize.x;
                else childNodeMin.x -= childNodeSize.x;
                if (nodeIndex & 2) childNodeMax.y += childNodeSize.y;
                else childNodeMin.y -= childNodeSize.y;
                if (nodeIndex & 4) childNodeMax.z += childNodeSize.z;
                else childNodeMin.z -= childNodeSize.z;

                node->children[nodeIndex] = std::make_shared<OctreeNode>(childNodeMin, childNodeMax);
            }

            InsertChunkIntoNode(chunk, node->children[nodeIndex], depth++);
        }
    }

    void BuildBoundingVolumes(std::shared_ptr<OctreeNode> node)
    {
        if (node->isLeaf)
        {
            for (uint32_t i = 0; i < node->chunks.size(); ++i)
            {
                node->bv.ExtendBy(node->chunks[i]->bv);
            }
        }

        else
        {
            for (uint8_t i = 0; i < 8; ++i)
            {
                if (node->children[i])
                {
                    BuildBoundingVolumes(node->children[i]);
                    node->bv.ExtendBy(node->children[i]->bv);
                }
            }
        }
    }

    uint8_t maxDepth = 6;
    std::shared_ptr<OctreeNode> root;
    const glm::vec3 sceneMin;
    const glm::vec3 sceneMax;
};

class BVH
{
public:
    BVH(std::vector<std::shared_ptr<TriMesh>> meshes) : meshes(meshes)
    {
        // Calculate scene extents and number of triangles
        uint32_t numTriangles = 0;
        glm::vec3 sceneMax(-Infinity);
        glm::vec3 sceneMin(Infinity);
        for (std::shared_ptr<TriMesh> mesh : meshes)
        {
            numTriangles += mesh->triangles.size();
            for (Triangle triangle : mesh->triangles)
            {
                sceneMax = glm::max(triangle.v0, sceneMax);
                sceneMin = glm::min(triangle.v0, sceneMin);

                sceneMax = glm::max(triangle.v1, sceneMax);
                sceneMin = glm::min(triangle.v1, sceneMin);

                sceneMax = glm::max(triangle.v2, sceneMax);
                sceneMin = glm::min(triangle.v2, sceneMin);
            }
        }

        glm::vec3 sceneSize = sceneMax - sceneMin;

        // Calculate grid resolution
        // Cleary et al. 1983. Design and analysis of a parallel ray tracing computer.
        // In their formula, lambda controls the granularity of the grid;
        // values between 3 and 5 are recommended
        float lambda = 3.f;
        float sceneVolume = sceneSize.x * sceneSize.y * sceneSize.z;
        // Calculate resolution of grid, and create a chunk per grid cell
        glm::vec3 resolution = glm::floor(sceneSize * glm::pow((glm::vec3(static_cast<float>(numTriangles)) / sceneVolume) * lambda, glm::vec3(1.f / 3.f)));

        resolution = glm::clamp(resolution, glm::vec3(1.f), glm::vec3(128.f));
        uint32_t numChunks = static_cast<uint32_t>(resolution.x * resolution.y * resolution.z);
        chunks = std::vector<std::shared_ptr<Chunk>>(numChunks + 1);

        // Add triangles to chunks
        for (uint32_t j = 0; j < meshes.size(); ++j)
        {
            for (uint32_t i = 0; i < meshes[j]->triangles.size(); ++i)
            {
                glm::vec3 triangleMin = glm::min(glm::min(meshes[j]->triangles[i].v0, meshes[j]->triangles[i].v1), meshes[j]->triangles[i].v2) - sceneMin;
                glm::vec3 chunkCoords = glm::floor((triangleMin / sceneSize) * resolution);
                uint32_t chunkIndex = static_cast<uint32_t>(glm::floor(chunkCoords.x * resolution.y * resolution.z + chunkCoords.y * resolution.z + chunkCoords.x));
                // TODO: Find source of overshoots
                chunkIndex = glm::min(chunkIndex, numChunks);
                if (!chunks[chunkIndex])
                {
                    std::shared_ptr<Chunk> chunk = std::make_shared<Chunk>();
                    chunks[chunkIndex] = chunk;
                }
                chunks[chunkIndex]->Insert(meshes[j], i);
            }
        }

        // Create octree structure from chunks
        octree = std::make_unique<Octree>(sceneMin, sceneMax);

        for (uint32_t i = 0; i < chunks.size(); ++i)
        {
            if (chunks[i])
            {
                octree->Insert(chunks[i]);
            }
        }

        // Create bounding volume for each octree node, from leaves to root
        octree->Build();
    }

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        return octree->Intersect(ray, isect);
    }

    std::vector<std::shared_ptr<Chunk>> chunks;
    std::vector<std::shared_ptr<TriMesh>> meshes;
    std::shared_ptr<Octree> octree;
};

std::shared_ptr<TriMesh> LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material)
{
    std::ifstream file;
    uint32_t numFaces;
    file.open(filePath);

    if (!file)
    {
        std::cerr << "Error: Mesh file " << filePath << " could not be opened.\n";
        return nullptr;
    }

    file >> numFaces;

    std::vector<uint32_t> faces(numFaces);
    uint32_t maxNormIndex = 0;

    uint32_t numVertIndices = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        if (!(file >> faces[i]))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }
        numVertIndices += faces[i];
    }

    std::vector<uint32_t> vertIndices(numVertIndices);
    uint32_t maxVertIndex = 0;

    uint32_t k = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i]; ++j)
        {
            uint32_t vertIndex;

            if (!(file >> vertIndex))
            {
                std::cerr << "Error: Mesh file could not be read.\n";
                return nullptr;
            }

            vertIndices[k] = vertIndex;

            k += 1;

            maxVertIndex = glm::max(maxVertIndex, vertIndex);
        }
    }

    uint32_t numVertCoords = (maxVertIndex + 1) * 3;
    std::vector<float> vertCoords(numVertCoords);

    for (uint32_t i = 0; i < numVertCoords; ++i)
    {
        float vertCoord;
        if (!(file >> vertCoord))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        vertCoords[i] = vertCoord;
    }

    std::vector<uint32_t> normIndices(numVertIndices);
    k = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i]; ++j)
        {
            uint32_t normIndex;

            if (!(file >> normIndex))
            {
                std::cerr << "Error: Mesh file could not be read.\n";
                return nullptr;
            }

            normIndices[k] = normIndex;
            k += 1;

            maxNormIndex = glm::max(maxNormIndex, normIndex);
        }
    }

    uint32_t numNormCoords = (maxNormIndex + 1) * 3;
    std::vector<float> normCoords(numNormCoords);

    for (uint32_t i = 0; i < numNormCoords; ++i)
    {
        float normCoord;
        if (!(file >> normCoord))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        normCoords[i] = normCoord;
    }



    // Calculate number of tris
    uint32_t numTris = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        numTris += faces[i] - 2;
    }

    // Create vector of vertices
    std::vector<glm::vec3> verts(maxVertIndex + 1);
    uint32_t j = 0;
    for (uint32_t i = 0; i < maxVertIndex + 1; ++i)
    {
        verts[i] = glm::transpose(objectToWorld) * glm::vec4(vertCoords[j], vertCoords[j + 1], vertCoords[j + 2], 1.f);
        j += 3;
    }

    // Create vector of normals
    std::vector<glm::vec3> norms(maxNormIndex + 1);
    j = 0;
    for (uint32_t i = 0; i < maxNormIndex + 1; ++i)
    {
        // Need to transform normals by inverse transpose - otherwise the normal will get scaled by any non-uniform scaling
        norms[i] = glm::normalize(glm::vec3(glm::inverse(((objectToWorld))) * glm::vec4(normCoords[j], normCoords[j + 1], normCoords[j + 2], 0.f)));
        j += 3;
    }

    // Create tri indices
    std::vector<uint32_t> triIndices(numTris * 3);
    k = 0;
    uint32_t l = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i] - 2; ++j)
        {
            triIndices[k] = vertIndices[l];
            triIndices[k + 1] = vertIndices[l + j + 1];
            triIndices[k + 2] = vertIndices[l + j + 2];
            k += 3;
        }
        l += faces[i];
    }

    // Create tri normal indices
    std::vector<uint32_t> triNormIndices(numTris * 3);
    k = 0;
    l = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i] - 2; ++j)
        {
            triNormIndices[k] = normIndices[l];
            triNormIndices[k + 1] = normIndices[l + j + 1];
            triNormIndices[k + 2] = normIndices[l + j + 2];
            k += 3;
        }
        l += faces[i];
    }

    std::shared_ptr<TriMesh> mesh = std::make_unique<TriMesh>(material);
    mesh->triangles.reserve(numTris);

    // Now we can finally build our triangles
    j = 0;
    for (uint32_t i = 0; i < numTris; ++i)
    {
        Triangle tri(verts[triIndices[j]], verts[triIndices[j + 1]], verts[triIndices[j + 2]], norms[triNormIndices[j]], norms[triNormIndices[j + 1]], norms[triNormIndices[j + 2]]);
        mesh->triangles.emplace_back(tri);

        j += 3;
    }

    return mesh;
}

struct Pixel
{
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

float StratifiedSample1D(std::default_random_engine& rng, uint32_t n, uint32_t nSamples)
{
    float invNSamples = 1.f / static_cast<float>(nSamples);
    std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());
    return (static_cast<float>(n) + distribution(rng)) * invNSamples;
}

// Using this technique instead of stratifying in x and y independently lets have any number of spp
// including primes, e.g 2, 3, 5, 7... and have them be evenly distributed
void LatinSquare(std::default_random_engine& rng, uint32_t nSamples, std::vector<glm::vec2>& samples)
{
    // Create stratified samples along diagonal
    for (uint32_t i = 0; i < nSamples; ++i)
    {
        samples[i] = glm::vec2(StratifiedSample1D(rng, i, nSamples), StratifiedSample1D(rng, i, nSamples));
    }

    // Shuffle dimensions
    for (uint32_t i = 0; i < nSamples; ++i)
    {
        std::uniform_int_distribution<uint32_t> distribution(0, nSamples - 1 - i);
        uint32_t choice = distribution(rng);
        std::swap(samples[i].x, samples[choice].x);
        choice = distribution(rng);
        std::swap(samples[i].y, samples[choice].y);
    }
}

float Gaussian(float width, float x)
{
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * glm::pi<float>() * sigma * sigma)) * glm::exp(-(x * x) / (2.f * sigma * sigma));
}

void AddSample(const RenderInfo& info, const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels)
{
    // Calculate discrete x and y bounds of filter
    uint32_t x0 = static_cast<uint32_t>(glm::floor(sampleCoords.x - info.filterWidth));
    uint32_t x1 = static_cast<uint32_t>(glm::ceil(sampleCoords.x + info.filterWidth));
    uint32_t y0 = static_cast<uint32_t>(glm::floor(sampleCoords.y - info.filterWidth));
    uint32_t y1 = static_cast<uint32_t>(glm::ceil(sampleCoords.y + info.filterWidth));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Calculate distance from sampleCoords
            float distX = (static_cast<float>(x) + 0.5f) - sampleCoords.x;
            float distY = (static_cast<float>(y) + 0.5f) - sampleCoords.y;
            float dist = glm::sqrt(distX * distX + distY * distY);

            // Get index into precomputed filter table based on distance from sample, as opposed to:
            // float filterWeight = Gaussian(scene.info.filterWidth, dist);
            uint8_t filterIndex = glm::min(static_cast<uint8_t>((dist / info.filterWidth) * FilterTableResolution), static_cast<uint8_t>(FilterTableResolution - 1));
            float filterWeight = filterTable[filterIndex];

            // Transform image coords into discrete tile coords
            uint32_t tileX = static_cast<uint32_t>(glm::floor(distX + glm::mod(sampleCoords.x - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));
            uint32_t tileY = static_cast<uint32_t>(glm::floor(distY + glm::mod(sampleCoords.y - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));

            uint32_t tileIndex = (tileY * info.tileSize) + tileX;

            pixels[tileIndex].contribution += L * filterWeight;
            pixels[tileIndex].filterWeightSum += filterWeight;
        }
    }
}

struct Camera
{
    Camera(float fov, glm::mat4 cameraToWorld) : fov(fov), cameraToWorld(cameraToWorld) {}

    const float fov = 10.f;
    const glm::mat4 cameraToWorld = glm::mat4(1.f);
};

class Scene
{
public:
    Scene(RenderInfo info, Camera camera, std::vector<std::shared_ptr<Light>> lights, std::shared_ptr<BVH> bvh) : info(info), camera(camera), lights(lights), bvh(bvh) {}

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        return bvh->Intersect(ray, isect);

        // for (std::shared_ptr<Light> light : lights)
        // {
        //     // Intersect each light
        // }
    }

    RenderInfo info;
    Camera camera;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<BVH> bvh;
};


Scene LoadScene(std::string scenePath)
{
    nlohmann::json json;
    std::ifstream ifs;
    ifs.open(scenePath);

    if (!ifs)
    {
        std::cerr << "Error: Scene file could not be opened.\n";
    }

    try
    {
        ifs >> json;
    }

    catch (nlohmann::json::exception& e)
    {
        std::cerr << "Error parsing scene: " << e.what() << "\nAborting.\n";
        std::abort();
    }

    nlohmann::json jsonInfo = json["renderInfo"];

    // Default values if not found in json
    uint32_t imageWidth = 64;
    uint32_t imageHeight = 64;
    uint32_t bucketSize = 32;
    uint32_t spp = 1;
    float filterWidth = 1.f;

    if (!jsonInfo.is_null())
    {
        try
        {
            imageWidth = jsonInfo["imageWidth"].get<uint32_t>();
            imageHeight = jsonInfo["imageHeight"].get<uint32_t>();
            bucketSize = jsonInfo["bucketSize"].get<uint32_t>();
            spp = jsonInfo["spp"].get<uint32_t>();
            filterWidth = jsonInfo["filterWidth"].get<float>();
        }

        catch (nlohmann::json::exception& e)
        {
            std::cerr << "Error in renderInfo: " << e.what() << "\n";
        }
    }

    else
    {
        std::cerr << "Warning: No render info found.\n";
    }

    RenderInfo info(imageWidth, imageHeight, bucketSize, spp, filterWidth);

    // Default values if not found in json
    glm::mat4 cameraToWorld(1.f);
    float fov = 1.f;
    std::vector<float> camTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };

    try {
        fov = json["camera"]["fov"].get<float>();
        camTransform = json["camera"]["transform"].get<std::vector<float>>();
    }

    catch (nlohmann::json::exception& e)
    {
        std::cerr << "Error in camera: " << e.what() << "\n";
    }

    uint8_t j = 0;
    for (uint8_t i = 0; i < 4; ++i)
    {
        cameraToWorld[i] = glm::vec4(camTransform[j], camTransform[j + 1], camTransform[j + 2], camTransform[j + 3]);
        j += 4;
    }

    Camera camera(fov, cameraToWorld);

    std::vector<std::string> mesh_filePaths;
    std::vector<std::shared_ptr<Material>> mesh_materials;
    if (!json["meshes"].is_null()) {
        for (auto& elem : json["meshes"])
        {
            std::string filePath = "";
            try
            {
                filePath = elem["filePath"].get<std::string>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh file path: " << e.what() << "\n";
            }
            mesh_filePaths.push_back(filePath);

            for (auto& mat : elem["material"])
            {
                std::shared_ptr<Material> material;
                try
                {
                    std::string type = mat["type"].get<std::string>();

                    if (type == "lambert")
                    {
                        std::vector<float> rhoGet = mat["rho"].get<std::vector<float>>();
                        glm::vec3 rho = glm::vec3(rhoGet[0], rhoGet[1], rhoGet[2]);
                        material = std::make_shared<DiffuseMaterial>(rho);
                    }

                    else if (type == "specular")
                    {
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        material = std::make_shared<SpecularMaterial>(R, eta);
                    }

                    else if (type == "glossy")
                    {
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        // float alpha = AlphaToRoughness(mat["roughness"].get<float>()); 
                        // float alpha = roughness * roughness;
                        // float alpha = RoughnessToAlpha(mat["roughness"].get<float>());
                        float alpha = mat["alpha"].get<float>();
                        material = std::make_shared<GlossyDielectricMaterial>(R, eta, alpha);
                    }

                    else if (type == "plastic")
                    {
                        std::vector<float> rhoGet = mat["rho"].get<std::vector<float>>();
                        glm::vec3 rho = glm::vec3(rhoGet[0], rhoGet[1], rhoGet[2]);
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        // float alpha = AlphaToRoughness(mat["roughness"].get<float>()); 
                        // float alpha = roughness * roughness;
                        // float alpha = RoughnessToAlpha(mat["roughness"].get<float>());
                        float alpha = mat["alpha"].get<float>();
                        material = std::make_shared<PlasticMaterial>(rho, R, eta, alpha);
                    }

                    else
                    {
                        std::cerr << "Error: '" << type << "' is not a material type.\nAborting.\n";
                        std::abort();
                    }
                }
                catch (nlohmann::json::exception& e)
                {
                    std::cerr << "Error in mesh material type: " << e.what() << "\n";
                }
                mesh_materials.push_back(material);
            }
        } 
    }

    else
    {
        std::cerr << "Error: No meshes found.\n";
    }

    std::vector<glm::mat4> meshTransforms;
    if (!json["meshes"].is_null()) {
        for (auto& elem : json["meshes"])
        {
            std::vector<float> meshTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };
            try
            {
                meshTransform = elem["transform"].get<std::vector<float>>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh transform: " << e.what() << "\n";
            }

            glm::mat4 objectToWorld(1.f);
            uint8_t j = 0;
            for (uint8_t i = 0; i < 4; ++i)
            {
                objectToWorld[i] = glm::vec4(meshTransform[j], meshTransform[j + 1], meshTransform[j + 2], meshTransform[j + 3]);
                j += 4;
            }
            meshTransforms.push_back(objectToWorld);
        }
    }

    std::vector<std::shared_ptr<TriMesh>> meshes;

    for (uint32_t i = 0; i < mesh_filePaths.size(); ++i) {
        meshes.push_back(LoadMeshFromFile(mesh_filePaths[i], meshTransforms[i], mesh_materials[i]));
    }

    std::cout << "Building BVH...\n";
    std::shared_ptr<BVH> bvh = std::make_shared<BVH>(BVH(meshes));

    std::vector<std::shared_ptr<Light>> lights;

    if (!json["lights"].is_null()) {
        for (auto& elem : json["lights"])
        {
            std::vector<float> lightTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };
            try
            {
                lightTransform = elem["transform"].get<std::vector<float>>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in light transform: " << e.what() << "\n";
            }

            glm::mat4 lightToWorld(1.f);
            uint8_t j = 0;
            for (uint8_t i = 0; i < 4; ++i)
            {
                lightToWorld[i] = glm::vec4(lightTransform[j], lightTransform[j + 1], lightTransform[j + 2], lightTransform[j + 3]);
                j += 4;
            }

            try
            {
                std::string type = elem["type"].get<std::string>();

                if (type == "distant")
                {
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<DistantLight>(Le, intensity, lightToWorld);
                    lights.push_back(light);
                }

                if (type == "disk")
                {
                    float radius = elem["radius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<DiskLight>(radius, Le, intensity, lightToWorld);
                    lights.push_back(light);
                }

                if (type == "ring")
                {
                    float radius = elem["radius"].get<float>();
                    float innerRadius = elem["innerRadius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<RingLight>(radius, innerRadius, Le, intensity, lightToWorld);
                    lights.push_back(light);
                }
            }

            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh material type: " << e.what() << "\n";
            }
        }
    }

    // Temp code
    // std::shared_ptr<Light> diskLight = std::make_shared<DiskLight>(glm::vec3(5.f, 5.f, 5.f), glm::mat4((1.0, 0.0, 0.0, 1.0000, 0.0, 1.0, 0.0, 2.0000, 0.0000, 0.0, 1.0, 3.0000, 0.0000, 0.0000, 0.0000, 1.000)));
    // lights.push_back(diskLight);

    return Scene(info, camera, lights, bvh);
}

std::vector<Pixel> RenderTile(const Scene& scene, const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1)
{
    std::vector<Pixel> pixels((scene.info.tileSize) * (scene.info.tileSize));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Don't render beyond total image extents when the current tile goes beyond them
            if (x < scene.info.totalWidth && y < scene.info.totalHeight)
            {
                // Create stratified image samples in a Latin square pattern
                std::default_random_engine rng;
                rng.seed(y * scene.info.totalWidth + x);
                std::vector<glm::vec2> imageSamples(scene.info.spp);
                LatinSquare(rng, scene.info.spp, imageSamples);

                for (uint32_t i = 0; i < scene.info.spp; ++i)
                {
                    glm::vec2 imageSample = imageSamples[i];

                    // Create ray in camera-space
                    float aspectRatio = static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.imageHeight);
                    float px = (((static_cast<float>(x) + imageSample.x) / static_cast<float>(scene.info.imageWidth)) * 2.f - 1.f) * glm::tan(glm::radians(scene.camera.fov)) * aspectRatio;
                    float py = (((static_cast<float>(y) + imageSample.y) / static_cast<float>(scene.info.imageHeight)) * -2.f + 1.f) * glm::tan(glm::radians(scene.camera.fov));

                    glm::vec4 o = glm::vec4(0.f, 0.f, 0.f, 1.f);
                    glm::vec4 d(px, py, -1.f, 0.f);
                    glm::normalize(d);

                    // Transform ray into world-space using camera transform
                    glm::mat4 camToWorld = scene.camera.cameraToWorld;

                    o = o * camToWorld;
                    d = d * camToWorld;

                    Ray ray(o, d);

                    // Radiance
                    glm::vec3 L(0.f, 0.f, 0.f);
                    float alpha = 0.f;
                    Intersection isect;
                    // Throughput
                    glm::vec3 beta(1.f);
                    BSDFFlags flags;

                    float roughnessOffset = 0.f;

                    for (uint32_t bounce = 0; bounce < 10; ++bounce)
                    {
                        // TODO: Move this inside of scene.Intersect()
                        float lightTMax = isect.tMax;
                        bool lightHit = false;
                        glm::vec3 Le(0.f);
                        for (const auto& light : scene.lights)
                        {
                            Intersection lightIsect;
                            glm::vec3 Li = light->Li(&lightIsect, ray.o, ray.d);
                            if (lightIsect.tMax < lightTMax)
                            {
                                Le = Li;
                                lightTMax = lightIsect.tMax;
                                isect.tMax = lightIsect.tMax;
                                lightHit = true;
                            }
                        }

                        if (bounce == 0)
                        {
                            if (lightHit == true)
                            {
                                L = Le;
                                break;
                            }
                        }

                        if (scene.Intersect(ray, &isect))
                        {
                            if (bounce == 0) alpha = 1.f;

                            std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());

                            BSDF bsdf = isect.material->CreateBSDF(isect.sn, roughnessOffset);
                            glm::vec3 wo = bsdf.ToLocal(-ray.d);

                            float shadowBias = 0.0001f;

                            glm::vec3 wi;
                            glm::vec3 Li;
                            float scatteringPdf = 0.f;
                            float lightingPdf = 0.f;

                            float numLights = static_cast<float>(scene.lights.size());
                            uint8_t lightIndex = static_cast<uint8_t>(glm::min(distribution(rng), 1.f - glm::epsilon<float>()) * numLights);
                            std::shared_ptr<Light> light = scene.lights[lightIndex];
                            glm::vec3 f(0.f);

                            // Compute direct light
#if BSDF_SAMPLING
                            scatteringPdf = 0.f;
                            glm::vec2 scatterSample(distribution(rng), distribution(rng));
                            f = bsdf.Sample_f(wo, &wi, scatterSample, &scatteringPdf, &flags);
                                
                            if (scatteringPdf > 0.f)
                            {
                                Ray shadowRay(isect.p + (isect.gn * shadowBias), bsdf.ToWorld(wi));
                                Intersection lightIsect;
                                Li = light->Li(&lightIsect, isect.p, bsdf.ToWorld(wi));
                                if (!scene.bvh->Intersect(shadowRay, &lightIsect))
                                {
                                    float weight = 1.f;
                                
                                    if (!(flags & SPECULAR))
                                    {
                                        // TODO: I should probably do this in one function
                                        lightingPdf = light->Pdf(&lightIsect, isect.p, bsdf.ToWorld(wi));
#if LIGHT_SAMPLING
                                        weight = (scatteringPdf * scatteringPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                                        if (lightingPdf > 0.f)
                                        {
                                            L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / scatteringPdf;
                                            L *= numLights;
                                        }
                                    }

                                    else
                                    {
                                        L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / scatteringPdf;
                                        L *= numLights;
                                    }
                                }
                            }
#endif
#if LIGHT_SAMPLING
                            glm::vec3 wiWorld;
                            lightingPdf = 0.f;
                            glm::vec2 lightSample(distribution(rng), distribution(rng));
                            // lightIsect.tMax used for checking shadows
                            Intersection lightIsect;
                            Li = light->Sample_Li(&lightIsect, isect.p, &wiWorld, lightSample, &lightingPdf);
                                
                            Ray shadowRay(isect.p + (isect.gn * shadowBias), wiWorld);
                            if (!scene.bvh->Intersect(shadowRay, &lightIsect) && lightingPdf > 0.f)
                            {
                                float weight = 1.f; 
                                wi = bsdf.ToLocal(wiWorld);
                                scatteringPdf = bsdf.Pdf(wo, wi);
                                if (scatteringPdf > 0.f)
                                {
                                    glm::vec3 f = bsdf.f(wo, wi);
#if BSDF_SAMPLING
                                    weight = (lightingPdf * lightingPdf) / (scatteringPdf * scatteringPdf + lightingPdf * lightingPdf);
#endif
                                    L += (f * Li * glm::max(wi.z, 0.f) * beta * weight) / lightingPdf;
                                    // L *= numLights;
                                }
                            }
#endif

                            // Spawn new ray
                            glm::vec2 scatteringSample(distribution(rng), distribution(rng));
                            scatteringPdf = 0.f;
                            f = bsdf.Sample_f(wo, &wi, scatteringSample, &scatteringPdf, &flags);
                            if ( scatteringPdf <= 0.f ) break;
                            if ( flags & DIFFUSE ) roughnessOffset += 0.5f;
                            beta *= (f / scatteringPdf) * glm::abs(wi.z);
                            // Transform to world
                            ray = Ray(isect.p + (isect.gn * shadowBias), bsdf.ToWorld(wi));

                            // Russian roulette
                            float q = glm::max((beta.x + beta.y + beta.z) * 0.33333f, 0.f);
                            
                            // std::cout << q << "\n";

                            if (bounce > 3)
                            {
                                if (q >= distribution(rng))
                                {
                                    beta /= q;
                                }

                                else break;
                            }

                            isect = Intersection();

                            // if (x == 360 && y == scene.info.imageHeight - 188)
                            // {
                            //     // std::cout << isect.p.x << " " << isect.p.y << " " << isect.p.z << "\n";
                            //     // for (uint8_t i = 0; i < 16; ++i)
                            //     // {
                            //     //     glm::vec3 temp = ray.o + (float)i * ray.d;
                            //     //     std::cout << temp.x << " " << temp.y << " " << temp.z << "\n";
                            //     // }
                            //     
                            //     // std::cout << wi.x << " " << wi.y << " " << wi.z << "\n";
                            //     L = glm::vec3(1.f, 0.f, 1.f);
                            // }
                        }

                        else break;
                    }

                    // Transform image sample to "total" image coords (image including filter bounds)
                    glm::vec2 sampleCoords = glm::vec2(static_cast<float>(x + scene.info.filterBounds) + imageSample.x, static_cast<float>(y + scene.info.filterBounds) + imageSample.y);
                    // Add sample to pixels within filter width
                    // Sample coords are respective to total image size
                    glm::vec4 Lalpha(L, alpha);
                    AddSample(scene.info, filterTable, sampleCoords, Lalpha, pixels);
                }
            }
        }
    }

    return pixels;
}

std::vector<Pixel> Render(const Scene& scene)
{
    std::vector<Pixel> pixels(scene.info.totalHeight * scene.info.totalWidth);

    uint32_t nBucketsX = uint32_t(glm::ceil(static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBucketsY = uint32_t(glm::ceil(static_cast<float>(scene.info.imageHeight) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBuckets = nBucketsX * nBucketsY;

    // Pre-compute filter values (only need to do this in 1D as filter is currently only isotropic)
    float filterTable[FilterTableResolution];
    for (uint8_t i = 0; i < FilterTableResolution; ++i)
    {
        filterTable[i] = Gaussian(FilterTableResolution - 1, i);
    }

    // Create tiles and render each one in parallel
    std::vector<Pixel> p(scene.info.tileSize * scene.info.tileSize);
    std::vector<std::vector<Pixel>> tiles(nBuckets, p);

    tbb::task_group tg;

    uint32_t nBucketsComplete = 0;

#if DEBUG_BUCKET
    for (uint32_t y = DEBUG_BUCKET_Y; y < DEBUG_BUCKET_Y + 1; ++y)
    {
        for (uint32_t x = DEBUG_BUCKET_X; x < DEBUG_BUCKET_X + 1; ++x)
        {
#else
    for (uint32_t y = 0; y < nBucketsY; ++y)
    {
        for (uint32_t x = 0; x < nBucketsX; ++x)
        {
#endif
            uint32_t index = y * nBucketsX + x;
            tg.run([&scene, &filterTable, &tiles, index, x, y, &nBucketsComplete, nBuckets]
                {
                    std::vector<Pixel> tile = RenderTile(scene, filterTable, scene.info.bucketSize * x, scene.info.bucketSize * (x + 1), scene.info.bucketSize * y, scene.info.bucketSize * (y + 1));
                    tiles[index] = tile;
                    nBucketsComplete++;
                    // TODO: Multi-threaded logging?
                    std::cout << "\r" << (int)glm::floor(((float)nBucketsComplete / (float)nBuckets) * 100.f) << "%   " << std::flush;
                });
        }
    }

    tg.wait();

    // Combine tiles into image
    std::cout << "\nCombining tiles into image...\n";

    for (uint32_t j = 0; j < nBucketsY; ++j)
    {
        for (uint32_t i = 0; i < nBucketsX; ++i)
        {
            for (uint32_t y = 0; y < scene.info.tileSize; ++y)
            {
                for (uint32_t x = 0; x < scene.info.tileSize; ++x)
                {
                    std::vector<Pixel> v = tiles[j * nBucketsX + i];
                    uint32_t pX = x + (i * scene.info.bucketSize);
                    uint32_t pY = y + (j * scene.info.bucketSize);
                    // Ignore tile pixels beyond image bounds
                    if (pX < (scene.info.imageWidth + scene.info.filterBounds) && pY < (scene.info.imageHeight + scene.info.filterBounds))
                    {
                        uint32_t pIndex = (pY * scene.info.totalWidth) + pX;
                        pixels[pIndex].contribution += v[(y * scene.info.tileSize) + x].contribution;
                        pixels[pIndex].filterWeightSum += v[(y * scene.info.tileSize) + x].filterWeightSum;
                    }
                }
            }
        }
    }

    return pixels;
}

void WriteImageToEXR(const RenderInfo& info, const std::vector<Pixel>& pixels, const char* filePath)
{
    // Write image to EXR
    Imf::Array2D<Imf::Rgba> imagePixels(info.imageHeight, info.imageWidth);

    // Ignore pixels beyond image bounds
    for (uint32_t y = info.filterBounds; y < (info.imageHeight + info.filterBounds); ++y)
    {
        for (uint32_t x = info.filterBounds; x < (info.imageWidth + info.filterBounds); ++x)
        {
            glm::vec4 result(0.f);

            result = pixels[y * info.totalWidth + x].contribution / pixels[y * info.totalWidth + x].filterWeightSum;

            imagePixels[y - info.filterBounds][x - info.filterBounds].r = result.r;
            imagePixels[y - info.filterBounds][x - info.filterBounds].g = result.g;
            imagePixels[y - info.filterBounds][x - info.filterBounds].b = result.b;
            imagePixels[y - info.filterBounds][x - info.filterBounds].a = result.a;
        }
    }

    Imf::RgbaOutputFile file(filePath, info.imageWidth, info.imageHeight, Imf::WRITE_RGBA);
    file.setFrameBuffer(&imagePixels[0][0], 1, info.imageWidth);
    file.writePixels(info.imageHeight);
}

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    
    std::cout << "Loading " << argv[1] << "...\n";
    Scene scene = LoadScene(argv[1]);
    
    std::cout << "Rendering...\n";
    std::vector<Pixel> image = Render(scene);
    
    std::cout << "Writing to " << argv[2] << "...\n";
    WriteImageToEXR(scene.info, image, argv[2]);
        
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << "Completed in " << duration.count() << "s\n";

    // std::default_random_engine rng;
    // rng.seed(1);
    // std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());
    // for (int i = 0; i < 1024; ++i)
    // {
    //     glm::vec2 sample(distribution(rng), distribution(rng));
    //     float pdf = 0.f;
    //     glm::vec3 res = CosineSampleHemisphere(sample, &pdf);
    //     std::cout << res.x << " " << res.y << " " << res.z << "\n";
    // }

    // std::default_random_engine rng;
    // rng.seed(0);
    // std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());
    // for (int i = 0; i < 8; ++i)
    // {
    //     glm::vec3 sampleA = glm::normalize(glm::vec3(distribution(rng), distribution(rng), distribution(rng)));
    //     glm::vec3 sampleB = glm::normalize(glm::vec3(distribution(rng), distribution(rng), distribution(rng)));
    // 
    //     BSDF bsdf(sampleA, 1);
    //     glm::vec3 v = sampleB;
    //     glm::vec3 vWorld = bsdf.ToWorld(v);
    //     glm::vec3 vLocal = bsdf.ToLocal(vWorld);
    //     std::cout << sampleA.x << " " << sampleA.y << " " << sampleA.z << "\n\n";
    //     std::cout << v.x << " " << v.y << " " << v.z << "\n";
    //     // std::cout << vWorld.x << " " << vWorld.y << " " << vWorld.z << "\n";
    //     std::cout << vLocal.x << " " << vLocal.y << " " << vLocal.z << "\n\n";
    // }

    // std::default_random_engine rng;
    // rng.seed(0);
    // std::uniform_real_distribution<float> distribution(0.f, 1.f - glm::epsilon<float>());
    // for (int i = 0; i < 512; ++i)
    // {
    //     glm::vec2 sample(distribution(rng), distribution(rng));
    //     float pdf;
    //     glm::vec2 p = UniformSampleRing(sample, &pdf);
    //     std::cout << p.x << " " << p.y << " " << 0.f << "\n";
    // }

    return 0;
}