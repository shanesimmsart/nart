#include "../../include/nart/patterns/texturepattern.h"

TexturePattern::TexturePattern(std::string filePath, bool isRoughness)
    : isRoughness(isRoughness) {
    Imf::RgbaInputFile file(filePath.c_str());
    Imath_3_1::Box2i dw = file.dataWindow();

    width = dw.max.x - dw.min.x + 1;
    height = dw.max.y - dw.min.y + 1;
    pixels.resizeErase(height, width);

    file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);
}

glm::vec3 TexturePattern::GetValue(const Intersection& isect) {
    float u = glm::min(glm::max(isect.st.x, 0.0001f), 0.9999f);
    float v = glm::min(glm::max(1.f - isect.st.y, 0.0001f), 0.9999f);
    int indexU = (int)((float)width * u);
    int indexV = (int)((float)height * v);
    float r = pixels[indexV][indexU].r;
    float g = pixels[indexV][indexU].g;
    float b = pixels[indexV][indexU].b;

    if (isRoughness) {
        r *= r;
        g *= g;
        b *= b;
    }

    return glm::vec3(r, g, b);
}