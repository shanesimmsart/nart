#pragma once

#include "../core/material.h"

class DiffuseMaterial : public Material {
 public:
  DiffuseMaterial(const glm::vec3& rho);

  BSDF CreateBSDF(const glm::vec3& n, float alphaTweak,
                  MemoryArena& memoryArena);

 private:
  // Eventually, material inputs will be replaced with patterns, that vary
  // depending on UVs, position, etc.
  const glm::vec3 rho;
};
