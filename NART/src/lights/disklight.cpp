#include "disklight.h"

DiskLight::DiskLight(float radius, glm::vec3 Le, float intensity, glm::mat4 LightToWorld) : Light(Le, intensity, LightToWorld), radius(radius)
{
    isDelta = false;
}

glm::vec3 DiskLight::Li(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) const
{
    if (Pdf(lightIsect, p, wi) > 0.f) return Le * intensity;

    else return glm::vec3(0.f);
}

glm::vec3 DiskLight::Sample_Li(Intersection* lightIsect, const glm::vec3& p, glm::vec3* wi, glm::vec2 sample, float* pdf) const
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

float DiskLight::Pdf(Intersection* lightIsect, const glm::vec3& p, const glm::vec3& wi) const
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