# nart
This project is where I put things I have learned about path tracing / rendering in general into practice.

Everything is written without referring to existing code or ray tracing resources, just what I can remember from what I've learned and have been able to derive by myself.

You can follow my progress [here](https://twitter.com/shaneasimms/status/1728665955797217295).

![glass](https://github.com/user-attachments/assets/8cbb833f-176d-4527-ba49-d5751709814e)
![veach](https://github.com/user-attachments/assets/34fafd04-04a3-483d-85e6-2a1732a8b190)

## Usage:
*nart* takes a JSON scene file and an output file name as input, and renders the given scene as an EXR image.

For example on Linux; running:

`$ ./build/nart input/scenes/glassSphere.json output/glassSphere`

will render the included glass sphere JSON scene and render it as `output/glassSphere_0.exr` (if there are multiple render sessions given in the scene file, each one will be suffixed 0, 1, 2, etc.).

There are also optional flags that can be used to override the various rendering parameters:
```-imageWidth / -i
-imageHeight / -i
-bucketSize / -b
-spp / -s
-filterWidth / -f
-rougheningFactor / -r
```

For example, running:

`$ ./build/nart input/scenes/glassSphere.json output/glassSphere -s 64`

will override the spp set in the JSON scene file to 64.

## Features:
  - JSON scene file description
  - Tiled, multi-threaded rendering
  - BVH acceleration structure
  - Gaussian image filtering
  - EXR file format
  - Lambertian, microfacet specular and microfacet transmission BRDFs
  - Analytic disk and ring area lights
  - Multiple importance sampling
  - Russian roulette
  - Roughening over paths to improve sampling of caustics

## Requirements:
Supported operating systems:
  - Linux
  - Windows

Tested compilers:
  - Clang 14
  - MSVC 19

## Dependencies:
  - C++ 14
  - [OpenEXR](https://openexr.com/en/latest/install.html#install)
  - [GLM](https://github.com/g-truc/glm)
  - [Intel TBB](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#onetbb)
  - [nlohmann/json](https://github.com/nlohmann/json)
