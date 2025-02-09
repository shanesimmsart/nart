#include "../../include/nart/lights/environmentlight.h"

EnvironmentLight::EnvironmentLight(PatternPtr&& Le, float intensity,
                                   const glm::mat4& LightToWorld)
    : Le(std::move(Le)), intensity(intensity), LightToWorld(LightToWorld) {
    isDelta = false;
}

glm::vec3 EnvironmentLight::Li(Intersection& lightIsect, const glm::vec3& p,
                               const glm::vec3& wi, float* pdf) const {
    float theta = glm::acos(wi.z);
    float phi = glm::atan(wi.y, wi.x) + glm::pi<float>();
    // if (wi.x < 0.f) phi += glm::pi<float>();
    if (phi > glm::two_pi<float>()) phi -= glm::two_pi<float>();
    if (phi < 0.f) phi += glm::two_pi<float>();

    lightIsect.st = glm::vec2(1.f - (phi * glm::one_over_two_pi<float>()),
                              1.f - (theta * glm::one_over_pi<float>()));

    if (pdf) *pdf = glm::one_over_pi<float>() * 0.25f;

    lightIsect.tMax = 0x7f7ffffef;

    return Le->GetValue(lightIsect) * intensity;
}

glm::vec3 EnvironmentLight::Sample_Li(Intersection& lightIsect,
                                      const glm::vec3& p, glm::vec3& wi,
                                      glm::vec2 sample, float& pdf) const {
    // wi = glm::vec3(0.f);
    // pdf = 0.f;
    // lightIsect.tMax = 0x7f7ffffef;
    // return glm::vec3(0.f);

    float theta = glm::acos(1.f - (2.f * sample.x));
    float phi = sample.y * 2.f * glm::pi<float>();

    float x = glm::cos(phi) * glm::sin(theta);
    float y = glm::sin(phi) * glm::sin(theta);
    float z = glm::cos(theta);

    wi = glm::vec3(x, y, z);

    lightIsect.st = glm::vec2(theta * glm::one_over_pi<float>(),
                              phi * glm::one_over_two_pi<float>());

    pdf = glm::one_over_pi<float>() * 0.25f;

    lightIsect.tMax = 0x7f7ffffef;

    return Le->GetValue(lightIsect) * intensity;
}

float EnvironmentLight::Pdf(Intersection& lightIsect, const glm::vec3& p,
                            const glm::vec3& wi) const {
    return 0.f;
}
