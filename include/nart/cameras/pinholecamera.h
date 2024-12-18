#pragma once

#include "../core/camera.h"

class PinholeCamera : public Camera
{
public:
    PinholeCamera(float fov, glm::mat4 cameraToWorld);
    Ray CastRay(const glm::vec2& imageSample, const uint32_t imageWidth, const uint32_t imageHeight, const uint32_t imageX, const uint32_t imageY) const;

private:
    const float fov = 10.f;
    const glm::mat4 cameraToWorld = glm::mat4(1.f);
};


