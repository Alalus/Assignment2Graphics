#include "SpotLight.h"

SpotLight::SpotLight(vec3 lightDirection) {
    this->lightType = Spot;
    this->lightDirection = lightDirection;
}