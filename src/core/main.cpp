#include <chrono>

#include "../../include/nart/core/nart.h"

int main(int argc, char* argv[]) {
    /*
    Imf::RgbaInputFile file("C:/Users/shane/Dev/nart/input/textures/uv.exr");
    Imath_3_1::Box2i dw = file.dataWindow();

    Imf::Array2D<Imf::Rgba> pixels;
    int width = dw.max.x - dw.min.x + 1;
    int height = dw.max.y - dw.min.y + 1;
    pixels.resizeErase(height, width);

    file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);

    for (uint8_t i = 0; i < 10; ++i)
    {
        float u = (float)i / 10.f;
        float v = 0.9f;
        int indexU = (int)((float)width * u);
        int indexV = (int)((float)height * v);
        // auto r = pixels[indexV][indexU].r;
        float r = pixels[indexV][indexU].r;
        float g = pixels[indexV][indexU].g;
        int n = 0;
    }
    */

    /*
    TexturePattern tp("C:/Users/shane/Dev/nart/input/textures/uv.exr");
    Intersection is;
    for (uint8_t i = 0; i < 10; ++i) {
        is.st = glm::vec2((float)i / 10.f, (float)i / 10.f);
        std::cout << tp.GetValue(is).y << std::endl;
    }
    */

    if (argc < 3) {
        std::cerr << "Too few arguments given.\nUsage example: " << argv[0]
                  << " <scene file> <output path>\n";
        return EXIT_FAILURE;
    }

    else {
        // Default render parameters
        // (Overridden by scene file if left at default)
        RenderParams params;
        if (!ParseRenderParamArguments(argc, argv, params)) {
            return EXIT_FAILURE;
        }

        std::cout << "Loading " << argv[1] << "...\n";
        Scene scene(argv[1]);
        std::vector<RenderSessionPtr> sessions =
            LoadSessions(argv[1], scene, params);

        if (sessions.empty()) {
            std::cerr << "Failed to load sessions from " << argv[1] << "\n";
            return EXIT_FAILURE;
        }

        uint8_t n = 0;
        for (RenderSessionPtr& session : sessions) {
            auto start = std::chrono::high_resolution_clock::now();

            std::cout << "Rendering...\n";
            std::vector<Pixel> image = session->Render();

            std::ostringstream oss;
            if (sessions.size() == 1)
                oss << argv[2] << ".exr";
            else
                oss << argv[2] << "_" << std::to_string(n++) << ".exr";
            std::string filePath = oss.str();

            std::cout << "Writing to " << filePath.c_str() << "...\n";
            session->WriteImageToEXR(image, filePath.c_str());

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> duration = end - start;
            std::cout << "Completed in " << duration.count() << "s\n";
        }
    }

    return 0;
}
