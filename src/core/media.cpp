#include "../../include/nart/core/media.h"

float DensityGrid::LookUp(uint8_t x, uint8_t y, uint8_t z) const {
    uint32_t index = (resolutionX * resolutionY * z) + (resolutionX * y) + x;

    return density[index];
}

// p is within [0, 1)^3
float DensityGrid::LookUp(glm::vec3 p) const {
    float x = glm::min(glm::max(0.f, p.x), 0.999f) *
              static_cast<float>(resolutionX - 1);
    uint8_t loX = x;
    uint8_t hiX = loX + 1;
    // Interpolation value
    float xD = (x - static_cast<float>(loX));

    float y = glm::min(glm::max(0.f, p.y), 0.999f) *
              static_cast<float>(resolutionY - 1);
    uint8_t loY = y;
    uint8_t hiY = loY + 1;
    // Interpolation value
    float yD = (y - static_cast<float>(loY));

    float z = glm::min(glm::max(0.f, p.z), 0.999f) *
              static_cast<float>(resolutionZ - 1);
    uint8_t loZ = z;
    uint8_t hiZ = loZ + 1;
    // Interpolation value
    float zD = (z - static_cast<float>(loZ));

    // Trilinearly interpolate value
    // 4 x lerps (xl to xh by xD)
    float x0 = glm::mix(LookUp(loX, loY, loZ), LookUp(hiX, loY, loZ), xD);
    float x1 = glm::mix(LookUp(loX, loY, hiZ), LookUp(hiX, loY, hiZ), xD);
    float x2 = glm::mix(LookUp(loX, hiY, loZ), LookUp(hiX, hiY, loZ), xD);
    float x3 = glm::mix(LookUp(loX, hiY, hiZ), LookUp(hiX, hiY, hiZ), xD);
    // 2 y lerps (yD between 4 x lerp outputs; 2 ly and 2 hy)
    float y0 = glm::mix(x0, x2, yD);
    float y1 = glm::mix(x1, x3, yD);
    // 1 z lerp  (zD between 2 y lerp outputs; 1 lz and 1 hz)
    float z0 = glm::mix(y0, y1, zD);

    return z0;
}

MajorantGrid::MajorantGrid(const DensityGrid& densityGrid, float sigma_maj)
    : sigma_maj(sigma_maj) {
    // Find maximum value out of density vertices + sampled majorant
    // vertices within density grid

    // Distance between majorant axes from 0 to densityGridRes-1
    float mx = ((static_cast<float>(densityGrid.resolutionX) - 1.f) / width);
    float my = ((static_cast<float>(densityGrid.resolutionY) - 1.f) / width);
    float mz = ((static_cast<float>(densityGrid.resolutionZ) - 1.f) / width);
    // Incremented distance between majorant axes
    float xi = mx;
    float yi = my;
    float zi = mz;

    uint8_t x0 = 0;
    uint8_t y0 = 0;
    uint8_t z0 = 0;
    uint8_t x1 =
        glm::min(static_cast<uint8_t>(glm::ceil(xi)), densityGrid.resolutionX);
    uint8_t y1 =
        glm::min(static_cast<uint8_t>(glm::ceil(yi)), densityGrid.resolutionY);
    uint8_t z1 =
        glm::min(static_cast<uint8_t>(glm::ceil(zi)), densityGrid.resolutionZ);

    for (uint8_t w = 0; w < width; ++w) {
        for (uint8_t v = 0; v < width; ++v) {
            for (uint8_t u = 0; u < width; ++u) {
                // Search through densityGrid vertex values within current
                // MajorantGrid voxel bounds
                float majorant = 0.f;
                for (uint8_t k = z0; k < z1; ++k) {
                    for (uint8_t j = y0; j < y1; ++j) {
                        for (uint8_t i = x0; i < x1; ++i) {
                            majorant =
                                glm::max(majorant, densityGrid.LookUp(i, j, k));
                        }
                    }
                }

                // Search through densityGrid values at MajorantGrid
                // vertices
                float px0 = static_cast<float>(u) * invWidth;
                float px1 = static_cast<float>(u + 1) * invWidth;
                float py0 = static_cast<float>(v) * invWidth;
                float py1 = static_cast<float>(v + 1) * invWidth;
                float pz0 = static_cast<float>(w) * invWidth;
                float pz1 = static_cast<float>(w + 1) * invWidth;

                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px0, py0, pz0)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px0, py0, pz1)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px0, py1, pz0)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px0, py1, pz1)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px1, py0, pz0)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px1, py0, pz1)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px1, py1, pz0)));
                majorant = glm::max(
                    majorant, densityGrid.LookUp(glm::vec3(px1, py1, pz1)));

                uint32_t index = (width * width * w) + (width * v) + u;
                majorants[index] = majorant * sigma_maj;

                xi += mx;
                x0 = x1;
                x1 = glm::min(static_cast<uint8_t>(glm::ceil(xi + 0.00001f)),
                              densityGrid.resolutionX);
            }
            yi += my;
            y0 = y1;
            y1 = glm::min(static_cast<uint8_t>(glm::ceil(yi + 0.00001f)),
                          densityGrid.resolutionY);
        }
        zi += mz;
        z0 = z1;
        z1 = glm::min(static_cast<uint8_t>(glm::ceil(zi + 0.00001f)),
                      densityGrid.resolutionZ);
    }
}

float MajorantGrid::LookUp(uint8_t x, uint8_t y, uint8_t z) const {
    uint32_t index = (width * width * z) + (width * y) + x;

    return majorants[index];
}

RayMajorantIterator::RayMajorantIterator(const Ray& ray,
                                         const glm::vec3& boundsMin,
                                         const glm::vec3& boundsMax,
                                         MajorantGrid* majorantGrid, float tMin,
                                         float tMax)
    : majorantGrid(majorantGrid),
      tMin(tMin),
      tCurrent(glm::max(0.f, tMin)),
      tMax(tMax) {
    // Calculate entering and exiting majorant voxel indices
    glm::vec3 pEntering = ray.o + (ray.d * tCurrent);
    glm::vec3 pExiting = ray.o + (ray.d * tMax);

    pEntering -= boundsMin;
    glm::vec3 boundsSize = boundsMax - boundsMin;
    pEntering /= boundsSize;
    pEntering.x = glm::max(glm::min(pEntering.x, 0.999999f), 0.f);
    pEntering.y = glm::max(glm::min(pEntering.y, 0.999999f), 0.f);
    pEntering.z = glm::max(glm::min(pEntering.z, 0.999999f), 0.f);
    pEntering *= majorantGrid->width;
    glm::ivec3 iEnter(pEntering);
    currentIndex = (majorantGrid->width * majorantGrid->width * iEnter.z) +
                   (majorantGrid->width * iEnter.y) + iEnter.x;

    pExiting -= boundsMin;
    pExiting /= boundsSize;
    pExiting.x = glm::max(glm::min(pExiting.x, 0.999999f), 0.f);
    pExiting.y = glm::max(glm::min(pExiting.y, 0.999999f), 0.f);
    pExiting.z = glm::max(glm::min(pExiting.z, 0.999999f), 0.f);
    pExiting *= majorantGrid->width;
    glm::ivec3 iExit(pExiting);
    finalIndex = (majorantGrid->width * majorantGrid->width * iExit.z) +
                 (majorantGrid->width * iExit.y) + iExit.x;

    // Calculate distance between axes along ray
    // Ray direction within majorant grids local [0-1)^3 space
    glm::vec3 gD = glm::normalize(pExiting - pEntering);
    if (pExiting == pEntering) gD = glm::vec3(1.f, 0.f, 0.f);
    crossDistance =
        glm::abs(((1.f / gD) * majorantGrid->invWidth) * boundsSize);
    if (gD.x == 0.f) crossDistance.x = Infinity;
    if (gD.y == 0.f) crossDistance.y = Infinity;
    if (gD.z == 0.f) crossDistance.z = Infinity;

    // Calculate next axis crossing distances along ray
    float tx = 0.f, ty = 0.f, tz = 0.f;
    if (ray.d.x >= 0.f) {
        tx = glm::abs((glm::ceil(pEntering.x + 0.00001f) - pEntering.x) / gD.x);
    } else {
        tx =
            glm::abs((glm::floor(pEntering.x - 0.00001f) - pEntering.x) / gD.x);
    }
    if (ray.d.y >= 0.f) {
        ty = glm::abs((glm::ceil(pEntering.y + 0.00001f) - pEntering.y) / gD.y);
    } else {
        ty =
            glm::abs((glm::floor(pEntering.y - 0.00001f) - pEntering.y) / gD.y);
    }
    if (ray.d.z >= 0.f) {
        tz = glm::abs((glm::ceil(pEntering.z + 0.00001f) - pEntering.z) / gD.z);
    } else {
        tz =
            glm::abs((glm::floor(pEntering.z - 0.00001f) - pEntering.z) / gD.z);
    }
    if (gD.x == 0.f) tx = Infinity;
    if (gD.y == 0.f) ty = Infinity;
    if (gD.z == 0.f) tz = Infinity;
    nextCrossing = glm::vec3(tx, ty, tz) * boundsSize * majorantGrid->invWidth;

    // Calculate whether each axis is positively or negatively stepped through
    stepAxis = glm::ivec3(1, 1, 1);
    if (gD.x < 0.f) stepAxis.x = -1;
    if (gD.y < 0.f) stepAxis.y = -1;
    if (gD.z < 0.f) stepAxis.z = -1;
}

bool RayMajorantIterator::Next(RayMajorantSegment& seg) {
    if (tCurrent + rayBias > tMax) {
        return false;
    }

    // Choose index of smallest next crossing
    // Truth table:
    // x < y    x < z   y < z   Index
    // ------------------------------
    // 0        0       0       2
    // 0        0       1       1
    // 0        1       0       0 (Not possible)
    // 0        1       1       1
    // 1        0       0       2
    // 1        0       1       0 (Not possible)
    // 1        1       0       0
    // 1        1       1       0
    uint8_t choice = 0;
    uint8_t choiceMap[8] = {2, 1, 0, 1, 2, 0, 0, 0};

    if (nextCrossing.x < nextCrossing.y) choice += 4;
    if (nextCrossing.x < nextCrossing.z) choice += 2;
    if (nextCrossing.y < nextCrossing.z) choice += 1;

    uint8_t index = choiceMap[choice];
    float dt = nextCrossing[index];

    seg = RayMajorantSegment(majorantGrid->majorants[currentIndex], tCurrent,
                             tCurrent + dt);

    if (done) return true;

    nextCrossing -= glm::vec3(dt);
    nextCrossing[index] = crossDistance[index];

    // TODO: This is gross
    currentIndex += glm::max(
        (int)glm::pow(majorantGrid->width, index) * stepAxis[index], 0);
    tCurrent += dt;

    return true;
}

float PhaseFunction::SamplePhaseFunction(glm::vec2 sample2D, glm::vec3& wi,
                                         float& pdf) {
    wi = UniformSampleSphere(sample2D, pdf);
    // Phase function is it's own PDF
    return pdf;
}

bool Medium::SampleMedium(glm::vec3(p), MediumProperties& mp) {
    if (p.x < boundsMin.x || p.y < boundsMin.y || p.z < boundsMin.z ||
        p.x > boundsMax.x || p.y > boundsMax.y || p.z > boundsMax.z) {
        return false;
    }

    // Scale p into [0-1)^3
    p -= boundsMin;
    p /= (boundsMax - boundsMin);

    float density = densityGrid.LookUp(p);

    mp = MediumProperties(PhaseFunction(), sigma_a * density, sigma_s * density,
                          Le * density);
    return true;
}

bool Medium::SampleRay(const Ray& ray, RayMajorantIterator& iter) {
    float tMin = -Infinity;
    float tMax = Infinity;
    float dotLow = Infinity;

    // Intersect AABB
    glm::vec3 normals[3];
    normals[0] = glm::vec3(1, 0, 0);
    normals[1] = glm::vec3(0, 1, 0);
    normals[2] = glm::vec3(0, 0, 1);

    for (uint8_t i = 0; i < 3; ++i) {
        // Equation of a plane:
        // (o+td)dotN - bound = 0
        // OdotN + tDdotN = bound
        // tDdotN = bound - OdotN
        // t = (bound - OdotN) / DdotN

        float slabMin = (boundsMin[i] - glm::dot(ray.o, normals[i])) /
                        glm::dot(ray.d, normals[i]);
        float slabMax = (boundsMax[i] - glm::dot(ray.o, normals[i])) /
                        glm::dot(ray.d, normals[i]);

        if (slabMin > slabMax) {
            std::swap(slabMin, slabMax);
        }

        if (slabMin > tMax || tMin > slabMax) {
            return false;
        }

        if (glm::dot(ray.o + (slabMin * ray.d), normals[i]) < dotLow) {
            dotLow = glm::dot(ray.o + (slabMin * ray.d), normals[i]);
        }

        tMin = glm::max(tMin, slabMin);
        tMax = glm::min(tMax, slabMax);
    }

    iter = RayMajorantIterator(ray, boundsMin, boundsMax, &majorantGrid, tMin,
                               tMax);

    return true;
}
