#pragma once

#include "reflection.h"

class Material
{
public:
    virtual BSDF CreateBSDF(glm::vec3 n, float roughnessOffset) = 0;

protected:
    Material() {};
};

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial(glm::vec3 rho);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    // Eventually, material inputs will be replaced with patterns, that vary depending on UVs, position, etc.
    glm::vec3 rho;
};

class SpecularMaterial : public Material
{
public:
    SpecularMaterial(glm::vec3 R, float eta);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    glm::vec3 R;
    float eta;
};

class GlossyDielectricMaterial : public Material
{
public:
    GlossyDielectricMaterial(glm::vec3 R, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    float alpha;
    glm::vec3 R;
    float eta;
};

class PlasticMaterial : public Material
{
public:
    PlasticMaterial(glm::vec3 rho, glm::vec3 R, float eta, float alpha);

    BSDF CreateBSDF(glm::vec3 n, float roughnessOffset);

private:
    glm::vec3 rho;
    float alpha;
    glm::vec3 R;
    float eta;
};


