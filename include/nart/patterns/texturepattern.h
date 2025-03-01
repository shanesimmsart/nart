#pragma once

#include "../core/pattern.h"
#include "../core/util.h"
// TODO: GET LOOKUP TABLE WORKING!!! >:d
//       OR F16C SSE
#define IMATH_HALF_USE_LOOKUP_TABLE 0
#define IMATH_HALF_NO_LOOKUP_TABLE
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfRgbaFile.h>

class Piecewise2DDistribution {
public:
    Piecewise2DDistribution(const Imf::Array2D<Imf::Rgba>& pixels,
                            uint32_t width,
          uint32_t height);

    glm::vec2 Sample(const glm::vec2& sample, float& pdf);

    float Pdf(const glm::vec2& sample);

private:
    uint32_t width = 0;
    uint32_t height = 0;
    float invW = 0.f;
    float invH = 0.f;
    std::vector<float> marginalPdf;
    std::vector<float> conditionalPdf;
    std::vector<float> marginalCdf;
    std::vector<float> conditionalCdf;
};

using Piecewise2DDistributionPtr = std::unique_ptr<Piecewise2DDistribution>;

class TexturePattern : public Pattern {
public:
    TexturePattern(std::string filePath, bool isRoughness = false, bool createPdf = false);

    glm::vec3 Sample(const glm::vec2& sample, glm::vec2& pdfSample, float& pdf);

    float Pdf(const glm::vec2& sample);

    glm::vec3 GetValue(const Intersection& isect);

private:
    uint16_t width, height;
    Imf::Array2D<Imf::Rgba> pixels;
    bool isRoughness = false;
    Piecewise2DDistributionPtr pdf2D = nullptr;
};
