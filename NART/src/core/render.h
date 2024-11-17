#pragma once

#include <tbb/task_group.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>

#include "scene.h"
#include "geometry.h"
#include "light.h"

struct Pixel
{
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

inline float Gaussian(float width, float x)
{
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * glm::pi<float>() * sigma * sigma)) * glm::exp(-(x * x) / (2.f * sigma * sigma));
}

class RenderSession
{
public:
    RenderSession(const Scene& scene, uint32_t imageWidth, uint32_t imageHeight, uint32_t bucketSize, uint32_t spp, float filterWidth);

    void AddSample(const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels);

    std::vector<Pixel> RenderTile(const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1);

    std::vector<Pixel> Render();

    void WriteImageToEXR(const std::vector<Pixel>& pixels, const char* filePath);

private:
    const Scene& scene;

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

std::vector<std::unique_ptr<RenderSession>> LoadSessions(std::string scenePath, const Scene& scene);


