#include <chrono>

#include "nart.h"

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    
    std::cout << "Loading " << argv[1] << "...\n";
    Scene scene = LoadScene(argv[1]);
    
    std::cout << "Rendering...\n";
    std::vector<Pixel> image = Render(scene);
    
    std::cout << "Writing to " << argv[2] << "...\n";
    WriteImageToEXR(scene.info, image, argv[2]);
        
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << "Completed in " << duration.count() << "s\n";

    return 0;
}


