
#pragma once
#include "glm/glm.hpp"
using namespace glm;

class Ray {
public:
    vec3 direction;
    vec3 position;

    Ray(vec3 direction, vec3 position) : direction(direction), position(position) {};

};

