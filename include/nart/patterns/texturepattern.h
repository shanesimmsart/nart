#pragma once

#include "../core/pattern.h"
// TODO: GET LOOKUP TABLE WORKING!!! >:d
//       OR F16C SSE
#define IMATH_HALF_USE_LOOKUP_TABLE 0
#define IMATH_HALF_NO_LOOKUP_TABLE
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfRgbaFile.h>

class TexturePattern : public Pattern {
public:
    TexturePattern(std::string filePath);

    glm::vec3 GetValue(const Intersection& isect);

private:
    uint16_t width, height;
    Imf::Array2D<Imf::Rgba> pixels;
};


