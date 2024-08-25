#pragma once

#include "geometry.h"

class Camera
{
public:
    virtual Ray CastRay(const glm::vec2& imageSample, const uint32_t& imageWidth, const uint32_t& imageHeight, const uint32_t& imageX, const uint32_t& imageY) const = 0;

    virtual ~Camera() {}

protected:
    Camera() {};
};



