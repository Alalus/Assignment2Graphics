#include "DirLight.h"
    
DirLight::DirLight(vec3 lightDirection){
    this->lightType = Directional;
    this->lightDirection = lightDirection;
}

