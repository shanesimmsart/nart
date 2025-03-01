#include "../../include/nart/patterns/texturepattern.h"

Piecewise2DDistribution::Piecewise2DDistribution(
    const Imf::Array2D<Imf::Rgba>& pixels, uint32_t width, uint32_t height)
    : width(width), height(height) {
    invW = 1.f / (float)width;
    invH = 1.f / (float)height;

    float fInt = 0.f;
    marginalPdf.resize(height);
    for (uint32_t j = 0; j < height; ++j) {
        marginalPdf[j] = 0.f;
        for (uint32_t i = 0; i < width; ++i) {
            marginalPdf[j] += std::abs(pixels[height - j - 1][i].r) +
                              std::abs(pixels[height - j - 1][i].g) +
                              std::abs(pixels[height - j - 1][i].b);
        }

        marginalPdf[j] *= invW;
        fInt += marginalPdf[j];
    }
    fInt /= (float)height;

    // Combining all 1D conditional PDFs for memory efficiency
    conditionalPdf.resize(width * height);
    for (uint32_t j = 0; j < height; ++j) {
        if (marginalPdf[j] != 0.f) {
            for (uint32_t i = 0; i < width; ++i) {
                conditionalPdf[j * width + i] =
                    std::abs(pixels[height - j - 1][i].r) +
                    std::abs(pixels[height - j - 1][i].g) +
                    std::abs(pixels[height - j - 1][i].b);
                conditionalPdf[j * width + i] /= marginalPdf[j];
            }
        }

        // Handle case where whole row of pixels is black
        else {
            for (uint32_t i = 0; i < width; ++i) {
                conditionalPdf[j * width + i] = 1.f;
            }
        }
    }

    float invFInt = 1.f / fInt;
    for (float& x : marginalPdf) {
        x *= invFInt;
    }

    marginalCdf.resize(height + 1);
    marginalCdf[0] = 0.f;
    marginalCdf[height] = 1.f;
    for (uint32_t i = 1; i < height; ++i) {
        marginalCdf[i] = marginalCdf[i - 1] + (marginalPdf[i - 1] * invH);
    }

    conditionalCdf.resize((width * height) + height);
    for (uint32_t i = 0; i < height; ++i) {
        conditionalCdf[i * (width + 1)] = 0.f;
        conditionalCdf[i * (width + 1) + width] = 1.f;
    }

    for (uint32_t j = 0; j < height; ++j) {
        for (uint32_t i = 1; i < width; ++i) {
            conditionalCdf[j * (width + 1) + i] =
                conditionalCdf[j * (width + 1) + i - 1] +
                (conditionalPdf[j * width + i - 1] * invW);
        }
    }
}

glm::vec2 Piecewise2DDistribution::Sample(const glm::vec2& sample, float& pdf) {
    // find CDF value below x using binary search
    uint32_t lowerBound =
        BinarySearch(sample.y, marginalCdf, 0, marginalCdf.size() - 1);

    float uc = 0.f;
    float vc =
        ((sample.y - marginalCdf[lowerBound]) / marginalPdf[lowerBound]) +
        ((float)lowerBound * invH);
    vc = glm::min(vc, 0.9999999f);

    uint32_t v = uint32_t(vc * height);
    if (marginalPdf[v] > 0.f) {
        lowerBound = BinarySearch(sample.x, conditionalCdf, v * (width + 1),
                                  v * (width + 1) + width);

        lowerBound %= (width + 1);
        uc = ((sample.x - conditionalCdf[v * (width + 1) + lowerBound]) /
              conditionalPdf[v * width + lowerBound]) +
             ((float)lowerBound * invW);
        uc = glm::min(uc, 0.9999999f);
        uint32_t u = uint32_t(uc * width);

        pdf = marginalPdf[v] * conditionalPdf[v * width + u];
    }

    else
        pdf == 0.f;

    return glm::vec2(uc, vc);
}

float Piecewise2DDistribution::Pdf(const glm::vec2& sample) {
    uint32_t u = uint32_t(sample.x * width);
    uint32_t v = uint32_t(sample.y * height);

    return marginalPdf[v] * conditionalPdf[v * width + u];
}

TexturePattern::TexturePattern(std::string filePath, bool isRoughness,
                               bool createPdf)
    : isRoughness(isRoughness) {
    Imf::RgbaInputFile file(filePath.c_str());
    Imath_3_1::Box2i dw = file.dataWindow();

    width = dw.max.x - dw.min.x + 1;
    height = dw.max.y - dw.min.y + 1;
    pixels.resizeErase(height, width);

    file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);

    if (createPdf) {
        pdf2D =
            std::make_unique<Piecewise2DDistribution>(pixels, width, height);
    }
}

glm::vec3 TexturePattern::Sample(const glm::vec2& sample, glm::vec2& pdfSample,
                                 float& pdf) {
    pdfSample = sample;

    if (!pdf2D) {
        pdf = 1.f;
    }

    else {
        pdfSample = pdf2D->Sample(sample, pdf);
    }

    float u = glm::min(glm::max(pdfSample.x, 0.0001f), 0.9999f);
    float v = glm::min(glm::max(1.f - pdfSample.y, 0.0001f), 0.9999f);
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

float TexturePattern::Pdf(const glm::vec2& sample) {
    if (!pdf2D) {
        return 1.f;
    }

    else
        return pdf2D->Pdf(sample);
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