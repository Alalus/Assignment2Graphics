#include "Plane.h"


Plane::Plane(vec4 details, objectTypes type) {
    this->objectDetails = details;
    this->objectType = type;
}

vec3 Plane::normal() {
    return vec3(
        objectDetails.x,
        objectDetails.y,
        objectDetails.z
    );
}

float Plane::distance() {
    return objectDetails.w;
}
float Plane::FindIntersection(Ray ray) {
    vec3 planeNormal = this->normal();

    // ax + by + cz + d = 0
    float a = planeNormal.x;
    float b = planeNormal.y;
    float c = planeNormal.z;
    float d = this->distance();

    float x0 = ray.position.x;
    float y0 = ray.position.y;
    float z0 = ray.position.z;

    float vx = ray.direction.x;
    float vy = ray.direction.y;
    float vz = ray.direction.z;

    // a*p_x(t) + b*p_y(t) + c*p_z(t) + d = 0
    // a(x0 + t * vx) + b(y0 + t * vy) + c(z0 + t * vz) + d = 0
    // => 
    float t = -(a * x0 + b * y0 + c * z0 + d) / (a * vx + b * vy + c * vz);
    return t;
}

vec3 Plane::getObjectColor(vec3 hitPoint) {
    // Checkerboard pattern
    float scale_parameter = 0.5f;
    float chessboard = 0;
    if (hitPoint.x < 0) {
        chessboard += floor((0.5 - hitPoint.x) / scale_parameter);
    }
    else {
        chessboard += floor(hitPoint.x / scale_parameter);
    }
    if (hitPoint.y < 0) {
        chessboard += floor((0.5 - hitPoint.y) / scale_parameter);
    }
    else {
        chessboard += floor(hitPoint.y / scale_parameter);
    }
    chessboard = (chessboard * 0.5) - int(chessboard * 0.5);
    chessboard *= 2;
    if (chessboard > 0.5) {
        return 0.5f * this->objectColorRGB;
    }
    return this->objectColorRGB;

}

vec3 Plane::getNormal(vec3 hitPoint) {
    return normalizedVector(normal());
}