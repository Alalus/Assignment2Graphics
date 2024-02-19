#pragma once

#include "glm/glm.hpp"
#include "Ray.h"
using namespace glm;

enum objectTypes {
    Regular,
    Transparent,
    Reflective,
    Space // backround
};
class Object {
public:
    objectTypes objectType;
    int objectIndex;
    vec4 objectDetails;
    vec3 objectColorRGB;
    float objectShiness;

    virtual void setObjectColor(vec4 color);
    virtual float FindIntersection(Ray ray) = 0;
    virtual vec3 getObjectColor(vec3 hitPoint) = 0;
    virtual vec3 getNormal(vec3 hitPoint) = 0;
};
vec3 normalizedVector(vec3 vec);

