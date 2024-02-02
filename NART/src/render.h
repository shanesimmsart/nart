#pragma once

#include <tbb/task_group.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>

#include "scene.h"
#include "geometry.h"
#include "lights.h"

struct Pixel
{
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

float Gaussian(float width, float x);

void AddSample(const RenderInfo& info, const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels);

std::vector<Pixel> RenderTile(const Scene& scene, const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1);

std::vector<Pixel> Render(const Scene& scene);

void WriteImageToEXR(const RenderInfo& info, const std::vector<Pixel>& pixels, const char* filePath);


