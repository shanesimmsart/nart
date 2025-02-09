#pragma once

#include <memory>
#include "memory.h"

struct Intersection;
class BSDF;

class Material {
public:
    // Creates BSDF and adds BxDFs
    // I think this can be const?
    // I've created a circular dependency!!! :d
    virtual BSDF CreateBSDF(const Intersection& isect, float alphaTweak,
                            MemoryArena& memoryArena) = 0;

    virtual ~Material() {}

protected:
    Material(){};
};

using MaterialPtr = std::unique_ptr<Material>;
