# nart
Not Another Ray Tracer is yet another one of my ray tracing projects.

Everything is written without referring to existing code or ray tracing resources, just what I can remember from what I've learned and have been able to derive by myself.

You can follow my progress [here](https://twitter.com/shaneasimms/status/1728665955797217295).

![glass](https://github.com/user-attachments/assets/8cbb833f-176d-4527-ba49-d5751709814e)
![progres 41](https://github.com/shanesimmsart/nart/assets/9335280/9f8b50ad-87e1-4efb-a1b8-a5a86f2ffd9a)

## Features so far:
  - JSON scene file description
  - Tiled, multi-threaded rendering
  - BVH acceleration structure
  - Gaussian image filtering
  - EXR file format
  - Lambertian and microfacet specular BRDFs
  - Analytic disk and ring area lights
  - Multiple importance sampling
  - Russian roulette
  - Firefly reduction via increasing roughness on secondary specular bounces after diffuse bounces (arguably better than clamping / filtering out caustics)

## Dependencies:
  - [OpenEXR](https://openexr.com/en/latest/install.html#install)
  - [GLM](https://github.com/g-truc/glm)
  - [Intel TBB](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#onetbb)
  - [nlohmann/json](https://github.com/nlohmann/json)
