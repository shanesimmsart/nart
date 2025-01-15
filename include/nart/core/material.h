#pragma once

#include "../bxdfs/dielectricbrdf.h"
#include "../bxdfs/lambertbrdf.h"
#include "../bxdfs/specularbrdf.h"
#include "../bxdfs/speculardielectricbrdf.h"
#include "../bxdfs/torrancesparrowbrdf.h"

class Material {
 public:
  // Creates BSDF and adds BxDFs
  // I think this can be const?
  virtual BSDF CreateBSDF(const glm::vec3& n, float alphaTweak,
                          MemoryArena& memoryArena) = 0;

  virtual ~Material() {}

 protected:
  Material(){};
};

using MaterialPtr = std::unique_ptr<Material>;
