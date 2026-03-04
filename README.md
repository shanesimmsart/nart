# nart
This project is where I put things that I have learned about path tracing into practice.

You can follow my progress [here](https://bsky.app/profile/shaneasimms.bsky.social/post/3lannomsuv327).

![lens1](https://github.com/user-attachments/assets/4d00b1e6-eb6e-41a6-be32-996d47974566)
<img width="480" height="854" alt="Nested dielectrics" src="https://github.com/user-attachments/assets/9cc25a22-1ccd-4918-b245-f04af29712f8" />
<img width="607" height="1079" alt="Volume integrator: scattering only" src="https://github.com/user-attachments/assets/5fdb5621-93f7-4c8a-8872-106f39e44437" />
<img width="607" height="1079" alt="Volume integrator: scattering, absorption and emission" src="https://github.com/user-attachments/assets/7819a95f-717c-40fd-92da-0b28ec37b086" />


## Usage:
**nart** takes a JSON scene file and an output file name as input, and renders the given scene as an EXR image.

For example on Linux; running:

`$ ./build/nart input/scenes/glassSphere.json output/glassSphere`


will render the included glass sphere JSON scene and write the result to `output/glassSphere.exr`.

There are also optional flags that can be used to override the various rendering parameters:
```
--imageWidth / -w
--imageHeight / -h
--bucketSize / -b
--spp / -s
--filterWidth / -f
--rougheningFactor / -r
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
  - Texture / normal mapping
  - Environment map importance sampling
  - Nested dielectrics

## Requirements:
Supported operating systems:
  - Linux
  - Windows

Tested compilers:
  - Clang 14
  - MSVC 19

## Dependencies:
  - C++14
  - [OpenEXR](https://openexr.com/en/latest/install.html#install)
  - [GLM](https://github.com/g-truc/glm)
  - [Intel TBB](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#onetbb)
  - [nlohmann/json](https://github.com/nlohmann/json)
