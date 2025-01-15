#include "../../include/nart/lights/distantlight.h"

DistantLight::DistantLight(const glm::vec3& Le, float intensity,
                           const glm::mat4& LightToWorld)
    : Le(Le), intensity(intensity), LightToWorld(LightToWorld) {
  isDelta = true;
  // Pointing down by default
  direction = glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld;
}

glm::vec3 DistantLight::Li(Intersection& lightIsect, const glm::vec3& p,
                           const glm::vec3& wi, float* pdf) const {
  if (pdf) *pdf = 0.f;
  return glm::vec3(0.f);
}

glm::vec3 DistantLight::Sample_Li(Intersection& lightIsect, const glm::vec3& p,
                                  glm::vec3& wi, glm::vec2 sample,
                                  float& pdf) const {
  wi = -direction;
  pdf = 1.f;
  return Le * intensity;
}

float DistantLight::Pdf(Intersection& lightIsect, const glm::vec3& p,
                        const glm::vec3& wi) const {
  return 0.f;
}
