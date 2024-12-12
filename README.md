# nart
Not Another Ray Tracer is yet another one of my ray tracing projects.

Everything is written without referring to existing code or ray tracing resources, just what I can remember from what I've learned and have been able to derive by myself.

You can follow my progress [here](https://twitter.com/shaneasimms/status/1728665955797217295).

![glass](https://github.com/user-attachments/assets/8cbb833f-176d-4527-ba49-d5751709814e)
![veach](https://github.com/user-attachments/assets/34fafd04-04a3-483d-85e6-2a1732a8b190)

## Features so far:
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

## Dependencies:
  - [OpenEXR](https://openexr.com/en/latest/install.html#install)
  - [GLM](https://github.com/g-truc/glm)
  - [Intel TBB](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#onetbb)
  - [nlohmann/json](https://github.com/nlohmann/json)
