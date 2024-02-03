#include <chrono>

#include "nart.h"

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    
    std::cout << "Loading " << argv[1] << "...\n";
    Scene scene = LoadScene(argv[1]);
    std::vector<std::shared_ptr<RenderSession>> sessions = LoadSessions(argv[1], scene);
    
    uint8_t n = 0;
    for (auto& session : sessions)
    {
        std::cout << "\nRendering...\n";
        std::vector<Pixel> image = session->Render();

        std::string seshID = std::to_string(n++);
        std::string extension = ".exr";
        std::string filePath = argv[2] + seshID + extension;

        std::cout << "Writing to " << filePath << "...\n";
        session->WriteImageToEXR(image, filePath.c_str());
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << "Completed in " << duration.count() << "s\n";

    return 0;
}


