#include "../../include/nart/core/geometry.h"

Ray::Ray(const glm::vec3& o, const glm::vec3& d) : o(o), d(d) {
    glm::vec3 absD = glm::abs(d);
    majorAxis = (absD.x > absD.y) ? ((absD.x > absD.z) ? 0 : 2)
                                  : ((absD.y > absD.z) ? 1 : 2);
    Sz = 1.f / d[majorAxis];
    uint8_t minorAxis0 = (majorAxis + 1);
    if (minorAxis0 == 3) minorAxis0 = 0;
    uint8_t minorAxis1 = (majorAxis + 2);
    if (minorAxis1 >= 3) minorAxis1 -= 3;
    Sx = -d[minorAxis0] * Sz;
    Sy = -d[minorAxis1] * Sz;
}

Triangle::Triangle(const glm::vec3& v0, const glm::vec3& v1,
                   const glm::vec3& v2, const glm::vec3& n0,
                   const glm::vec3& n1, const glm::vec3& n2,
                   const glm::vec2& uv0, const glm::vec2& uv1,
                   const glm::vec2& uv2)
    : v0(v0),
      v1(v1),
      v2(v2),
      n0(n0),
      n1(n1),
      n2(n2),
      uv0(uv0),
      uv1(uv1),
      uv2(uv2) {}

bool Triangle::Intersect(const Ray& ray, Intersection& isect) const {
    glm::vec3 n = glm::cross(v1 - v0, v2 - v0);
    // Back-face culling
    // if (glm::dot(n, ray.d) >= 0.f) return false;

    float t = (glm::dot(v0, n) - glm::dot(ray.o, n)) / glm::dot(ray.d, n);

    if (t <= isect.tMin || t >= isect.tMax) return false;

    // Translate ray origin to world origin
    glm::vec3 p0 = v0 - ray.o;
    glm::vec3 p1 = v1 - ray.o;
    glm::vec3 p2 = v2 - ray.o;

    // Permute axes so major axis of ray direction lies on +Z
    // glm::vec3 absD = glm::abs(ray.d);
    uint8_t majorAxis = ray.majorAxis;
    uint8_t minorAxis0 = (majorAxis + 1);
    if (minorAxis0 == 3) minorAxis0 = 0;
    uint8_t minorAxis1 = (majorAxis + 2);
    if (minorAxis1 >= 3) minorAxis1 -= 3;

    p0 = glm::vec3(p0[minorAxis0], p0[minorAxis1], p0[majorAxis]);
    p1 = glm::vec3(p1[minorAxis0], p1[minorAxis1], p1[majorAxis]);
    p2 = glm::vec3(p2[minorAxis0], p2[minorAxis1], p2[majorAxis]);

    // Shear ray direction to lie on +z
    // x' = x*1  + y*0  + z * (-dx/dz)
    // z' = x*0  + y*0  + z / dz
    // float Sz = 1.f / ray.d[majorAxis];
    p0.x += p0.z * ray.Sx;
    p0.y += p0.z * ray.Sy;
    p1.x += p1.z * ray.Sx;
    p1.y += p1.z * ray.Sy;
    p2.x += p2.z * ray.Sx;
    p2.y += p2.z * ray.Sy;

    // Check for intersection
    // Because p = (0, 0), cross products are simplified
    // u = distance from v1v2 (v0 is max)
    // v = distance from v2v0 (v1 is max)
    float e0 = (p1.x * p2.y) - (p1.y * p2.x);
    float e1 = (p2.x * p0.y) - (p2.y * p0.x);
    float e2 = (p0.x * p1.y) - (p0.y * p1.x);
    // We need to check if always on left OR right of edges, as permutation
    // can lead to inverse winding
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;

    if (glm::abs(e0) + glm::abs(e1) + glm::abs(e2) == 0.f) return false;

    // If intersection, normalise ray direction
    p0.z = p0.z * ray.Sz;
    p1.z = p1.z * ray.Sz;
    p2.z = p2.z * ray.Sz;

    // Calculate intersection
    float invDet = 1.f / (e0 + e1 + e2);
    glm::vec3 p = (v0 * e0 + v1 * e1 + v2 * e2) * invDet;

    float u = e0 * invDet;
    float v = e1 * invDet;

    n = glm::normalize(n);

    isect.tMax = t;
    isect.u = u;
    isect.v = v;
    isect.gn = n;
    isect.sn = n0 * isect.u + n1 * isect.v + n2 * (1 - isect.u - isect.v);
    isect.p = p;
    isect.st = uv0 * isect.u + uv1 * isect.v + uv2 * (1 - isect.u - isect.v);

    // TODO: Deal with 0 determinant
    float UVDet = (((uv0.x - uv2.x) * (uv1.y - uv2.y)) -
                   ((uv0.y - uv2.y) * (uv1.x - uv2.x)));
    float invUVDet = 1.f / UVDet;
    isect.dpds =
        ((v0 - v2) * (uv1.y - uv2.y) + (v1 - v2) * (uv2.y - uv0.y)) * invUVDet;
    isect.dpdt =
        ((v0 - v2) * (uv2.x - uv1.x) + (v1 - v2) * (uv0.x - uv2.x)) * invUVDet;

    return true;
}

TriMesh::TriMesh(std::vector<Triangle>&& tris, MaterialPtr&& material,
                 uint32_t meshID, uint8_t priority)
    : triangles(std::move(tris)), material(std::move(material)), meshID(meshID), priority(priority) {
    numTriangles = triangles.size();
}

const Triangle& TriMesh::GetTriangle(uint32_t index) const {
    if (index < numTriangles) {
        return triangles[index];
    }
}

const uint32_t TriMesh::GetNumTriangles() const { return numTriangles; }

const uint32_t TriMesh::GetMeshID() const { return meshID; }

const uint8_t TriMesh::GetPriority() const { return priority; }

Material* TriMesh::GetMaterial() const { return material.get(); }
