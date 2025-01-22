#include "../../include/nart/lights/disklight.h"

DiskLight::DiskLight(float radius, PatternPtr&& Le, float intensity,
                     const glm::mat4& LightToWorld)
    : radius(radius),
      Le(std::move(Le)),
      intensity(intensity),
      LightToWorld(LightToWorld) {
    isDelta = false;
}

glm::vec3 DiskLight::Li(Intersection& lightIsect, const glm::vec3& p,
                        const glm::vec3& wi, float* pdf) const {
    float LiPdf = Pdf(lightIsect, p, wi);

    if (LiPdf > 0.f) {
        if (pdf) *pdf = LiPdf;
        return Le->GetValue(lightIsect) * intensity;
    }

    else
        return glm::vec3(0.f);
}

glm::vec3 DiskLight::Sample_Li(Intersection& lightIsect, const glm::vec3& p,
                               glm::vec3& wi, glm::vec2 sample,
                               float& pdf) const {
    // Sample disk
    glm::vec4 diskSample =
        glm::vec4(UniformSampleDisk(sample) * radius, 0.f, 1.f);

    float u = ((diskSample.x + 1.f) * 0.5f) / radius;
    float v = ((diskSample.y + 1.f) * 0.5f) / radius;
    lightIsect.st = glm::vec2(u, 1.f - v);

    // Transform disk sample to world space
    diskSample = diskSample * LightToWorld;
    glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

    // Set w
    wi = glm::vec3(diskSample) - p;
    float distPToSample = glm::sqrt(wi.x * wi.x + wi.y * wi.y + wi.z * wi.z);
    wi = glm::normalize(wi);

    // Calculate pdf with respect to disk area
    pdf = 1.f / (glm::pi<float>() * radius * radius);

    float wiDotN = glm::dot(-wi, n);
    if (wiDotN <= 0.f) {
        pdf = 0.f;
        return glm::vec3(0.f);
    }
    // Calculate pdf projected onto hemisphere around p
    pdf = pdf * ((distPToSample * distPToSample) / wiDotN);

    lightIsect.p = diskSample;
    lightIsect.tMax = distPToSample;

    return Le->GetValue(lightIsect) * intensity;
}

float DiskLight::Pdf(Intersection& lightIsect, const glm::vec3& p,
                     const glm::vec3& wi) const {
    glm::vec3 diskCenter =
        glm::vec3(glm::vec4(0.f, 0.f, 0.f, 1.f) * LightToWorld);
    glm::vec3 n = glm::vec3(glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld);

    if (glm::dot(wi, n) >= 0.f) return 0.f;

    // Distance of plane from origin
    float D = glm::dot(diskCenter, n);

    // Intersect plane
    float t = (D - glm::dot(p, n)) / glm::dot(wi, n);
    if (t < 0.f) return 0.f;

    glm::vec3 pHit = p + t * wi;

    glm::vec3 diskCenterToPHit = pHit - diskCenter;

    float u = glm::dot(glm::vec4(diskCenterToPHit, 0.f),
                       glm::vec4(1.f, 0.f, 0.f, 0.f) * LightToWorld) /
              radius;
    float v = glm::dot(glm::vec4(diskCenterToPHit, 0.f),
                       glm::vec4(0.f, 1.f, 0.f, 0.f) * LightToWorld) /
              radius;
    u = (u + 1.f) * 0.5f;
    v = (v + 1.f) * 0.5f;
    lightIsect.st = glm::vec2(u, 1.f - v);

    // If not hit, pdf == 0
    float dist = diskCenterToPHit.x * diskCenterToPHit.x +
                 diskCenterToPHit.y * diskCenterToPHit.y +
                 diskCenterToPHit.z * diskCenterToPHit.z;
    if (dist > radius * radius) return 0.f;

    // If hit, calculate pdf projected onto hemisphere around p
    float pdf = 1.f / (glm::pi<float>() * radius * radius);
    pdf = pdf * ((t * t) / glm::dot(-wi, n));

    lightIsect.tMax = t;

    return pdf;
}
