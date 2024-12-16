#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "material.h"

#define Infinity std::numeric_limits<float>::infinity()

class Ray
{
public:
    Ray(glm::vec3 o, glm::vec3 d);

    glm::vec3 o = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 d = glm::vec3(0.f, 0.f, 1.f);

    uint8_t majorAxis = 2;
    float Sx = 1.f;
    float Sy = 1.f;
    float Sz = 1.f;
};

struct Intersection
{
    std::shared_ptr<Material> material;
    // Barycentric coords
    float u, v;
    // Point, geometric normal, shading normal
    glm::vec3 p, gn, sn;

    // Upper and lower bounds of ray to be considered
    float tMin = 0.f;
    float tMax = Infinity;

    // Flag for area light intersection
    bool isLight = false;
};

class Triangle
{
public:
    Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2);

    bool Intersect(const Ray& ray, Intersection* isect) const;

    const glm::vec3 v0, v1, v2;
    const glm::vec3 n0, n1, n2;
};

class TriMesh
{
public:
    TriMesh(std::shared_ptr<Material> material);

    // Triangles are passed directly, outside of the ctor, to reduce memory allocations
    std::vector<Triangle> triangles;
    std::shared_ptr<Material> material;
};

using TriMeshPtr = std::unique_ptr<TriMesh>;


