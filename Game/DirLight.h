#pragma once

#include "Light.h"

class DirLight : public Light {
public:
    DirLight(vec3 lightDirection);
};


