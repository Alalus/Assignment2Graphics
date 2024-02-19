#include "Light.h"

void Light::setIntensity(vec4 lightIntensity) {
    this->lightRgbIntensity = vec3(lightIntensity.r, lightIntensity.g, lightIntensity.b);
    this->lightShiness = lightIntensity.a;
}