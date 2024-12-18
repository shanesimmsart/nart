#include "../../include/nart/cameras/pinholecamera.h"

PinholeCamera::PinholeCamera(float fov, glm::mat4 cameraToWorld) : fov(fov), cameraToWorld(cameraToWorld) {}

Ray PinholeCamera::CastRay(const glm::vec2& imageSample, const uint32_t imageWidth, const uint32_t imageHeight, const uint32_t imageX, const uint32_t imageY) const
{
    // Create ray in camera-space
    float aspectRatio = static_cast<float>(imageWidth) / static_cast<float>(imageHeight);
    float px = (((static_cast<float>(imageX) + imageSample.x) / static_cast<float>(imageWidth)) * 2.f - 1.f) * glm::tan(glm::radians(fov)) * aspectRatio;
    float py = (((static_cast<float>(imageY) + imageSample.y) / static_cast<float>(imageHeight)) * -2.f + 1.f) * glm::tan(glm::radians(fov));

    glm::vec4 o = glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec4 d(px, py, -1.f, 0.f);
    glm::normalize(d);

    // Transform ray into world-space using camera transform
    o = o * cameraToWorld;
    d = d * cameraToWorld;

    return Ray(o, d);
}


