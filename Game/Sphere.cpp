#include "Sphere.h"

Sphere::Sphere(vec4 details, objectTypes type) {
    this->objectDetails = details;
    this->objectType = type;
}

vec3 Sphere::center() {
    return vec3(objectDetails.x, objectDetails.y, objectDetails.z);
}

float Sphere::radius() {
    return objectDetails.w;
}

float Sphere::FindIntersection(Ray ray) {
    glm::vec3 p0_c = ray.position - center(); // vector from sphere center c to p0 : p0 - c
    float a = glm::dot(ray.direction, ray.direction); // v^2
    float b = 2.0f * glm::dot(ray.direction, p0_c); // 2v * (p0 - c)
    float c = glm::dot(p0_c, p0_c) - radius() * radius(); // (p0 - c)^2 - r^2

    float discriminant = b * b - 4 * a * c;

    if (discriminant <= 0) {
        // no intersection
        return -1;
    }
    else {
        //  nearest t>0 value of the intersection, if exists
        float discriminant_sqrt = sqrt(discriminant);
        float t1 = (-b - discriminant_sqrt) / (2 * a);
        float t2 = (-b + discriminant_sqrt) / (2 * a);

        if (t1 > 0 && t2 > 0) {
            return glm::min(t1, t2); // The ray intersects the sphere, return the closest point
        }
        else if (t1 > 0) {
            return t1; // only t1 is valid
        }
        else if (t2 > 0) {
            return t2; // only t2 is valid
        }
        else {
            // The sphere is behind the ray or the ray starts inside the sphere
            return -1;
        }

    }
}

vec3 Sphere::getObjectColor(vec3 hitPoint) {
    return this->objectColorRGB;
}

vec3 Sphere::getNormal(vec3 hitPoint) {
    return normalizedVector(hitPoint - center());
}