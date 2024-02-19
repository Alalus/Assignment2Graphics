#pragma once

#include "Object.h"
#include "glm/glm.hpp"
using namespace glm;

class Hit {

public:
    vec3 hitPoint;
    Object* object;
    Hit(vec3 hitPoint, Object* object) : hitPoint(hitPoint), object(object) {};
    
};

