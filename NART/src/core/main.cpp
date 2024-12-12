#include <chrono>

#include "nart.h"

int main(int argc, char* argv[])
{ 
    if (argc < 3)
    {
        std::cerr << "Too few arguments given.\nUsage example: " << argv[0] << " <scene file> <output path>";
    }

    else
    {
        std::cout << "Loading " << argv[1] << "...\n";
        Scene scene(argv[1]);
        std::vector<RenderSessionPtr> sessions = LoadSessions(argv[1], scene);

        uint8_t n = 0;
        for (auto& session : sessions)
        {
            auto start = std::chrono::high_resolution_clock::now();

            std::cout << "Rendering...\n";
            std::vector<Pixel> image = session->Render();

            std::ostringstream oss;
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


