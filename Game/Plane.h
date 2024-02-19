#pragma once

#include "Object.h"

class Plane : public Object {

public:
    Plane(vec4 details, objectTypes type);
    vec3 normal();
    float distance();
    float FindIntersection(Ray ray);
    vec3 getObjectColor(vec3 hitPoint);
    vec3 getNormal(vec3 hitPoint);
};

