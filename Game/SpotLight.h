#pragma once
#include "Light.h"

class SpotLight : public Light {
public:
    SpotLight(vec3 lightDirection);
};

