#include "Config.h"
#include "Sphere.h"
#include "Plane.h"
#include "SpotLight.h"
#include "DirLight.h"
#include "Ray.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath> 


class SceneParser {
public:
    // parseSceneFile is a function that parses the scene file and returns the scene data
    static std::vector<std::vector<std::string>> parseSceneFile(const std::string& fileName) {
        std::ifstream sceneFile(fileName);//open the file
        std::string currentLine;//current line in the file
        std::vector<std::vector<std::string>> sceneData;//2d vector to store the scene data

        while (getline(sceneFile, currentLine)) {//read the file line by line till the end of the file
            std::stringstream textLineStream(currentLine);//convert the string to a stream
            std::string currentArgInLine;//current argument in the line
            std::vector<std::string> sceneLineArgs;//scene line arguments

            while (getline(textLineStream, currentArgInLine, ' ')) {//read the line by space and store the arguments in the sceneLineArgs
                sceneLineArgs.push_back(currentArgInLine);//store the arguments in the sceneLineArgs vector
            }
            if (!sceneLineArgs.empty()) {//if the sceneLineArgs is not empty then store the sceneLineArgs in the sceneData
                sceneData.push_back(sceneLineArgs);
            }
        }

        return sceneData;
    }
};
// readSceneFile is a function that reads the scene file and creates the scene 
void Config::readSceneFile(const std::string& fileName, int width, int height) {
    // parseSceneFile is a function that parses the scene file and returns the scene data
    auto sceneData = SceneParser::parseSceneFile(fileName);
    // processSceneData is a function that processes the scene data and sets the scene data accordingly like give the object a color, set the spotlight positions, etc.
    processSceneData(sceneData);
    // initializeScene is a function that initializes the scene with the given width and height
    initializeScene(width, height);
}


void Config::processSceneData(const std::vector<std::vector<std::string>>& sceneData) {
    for (const auto& lineData : sceneData) {
        if (lineData.size() < 5) continue; // Ensure lineData has at least 5 elements

        vec4 dataVector = vec4(std::stof(lineData[1]), std::stof(lineData[2]),
            std::stof(lineData[3]), std::stof(lineData[4]));

        processLineType(lineData[0], dataVector);
    }

    setColorForObjects();
    setSpotlightPositions();
}
// processLineType is a function that processes the line type and the data vector and sets the scene data accordingly 
void Config::processLineType(const std::string& lineType, const vec4& dataVector) {
    if (lineType == "e") { eye = vec3(dataVector.x, dataVector.y, dataVector.z); }
    else if (lineType == "a") { ambient = dataVector; }
    else if (lineType == "d") { addLight(dataVector); }
    else if (lineType == "p") { positions.push_back(dataVector); }
    else if (lineType == "i") { intensities.push_back(dataVector); }
    else if (lineType == "o") { addObject(dataVector, Regular); }
    else if (lineType == "r") { addObject(dataVector, Reflective); }
    else if (lineType == "t") { addObject(dataVector, Transparent); }
    else if (lineType == "c") { colors.push_back(dataVector); }
}
//e - camera position - (X,Y,Z) - camera position 4 value for tranfer additional data
//a - ambient light color - (R,G,B,A)
//o -regular object sphere or plane -(x,y,z) if "-" in the 4 value then it is a sphere, if not then it is a plane
//r - reflctive object sphere or plane 
//t - transparent object sphere or plane
//c - obejct color - (R,G,B,A) 
//d - light direction - (X,Y,Z) - if  the 4 value is 0.0 its directional light, if not then 1.0 it is a spot light
//p - light position - (X,Y,Z) - light position 
//i - light color - (R,G,B,A) - light color

void Config::addLight(const vec4& dataVector) {
    dataVector.w == 1.0
        ? lights.push_back(new SpotLight(vec3(dataVector.x, dataVector.y, dataVector.z)))
        : lights.push_back(new DirLight(vec3(dataVector.x, dataVector.y, dataVector.z)));
}

void Config::addObject(const vec4& dataVector, objectTypes type) {
    dataVector.w > 0 ? objects.push_back(new Sphere(dataVector, type))
        : objects.push_back(new Plane(dataVector, type));
}

void Config::setColorForObjects() {
    for (size_t i = 0; i < objects.size(); ++i) {
        objects[i]->objectIndex = i;
        objects[i]->setObjectColor(colors[i]);
    }
}

// setSpotlightPositions is a function that sets the spotlight positions
void Config::setSpotlightPositions() {
    size_t j = 0;
    for (size_t i = 0; i < lights.size(); ++i) {
        auto& light = lights[i];
        if (light->lightType == Spot && j < positions.size()) {
            vec3 point = vec3(positions[j].x, positions[j].y, positions[j].z);
            float cosAngle = positions[j].w;
            ++j;
            light->lightPosition = point;
            light->lightCosAngle = cosAngle;
        }
        light->setIntensity(intensities[i]);
    }
}
// initializeScene is a function that initializes the scene with the given width and height its map the pixel width and height to the image width and height between -1 and 1
void Config::initializeScene(int width, int height) {
    //it takes the width and height of the image and sets the pixel width and height accordingly

    pixelWidth = 2.0 / static_cast<float>(width);
    pixelHeight = 2.0 / static_cast<float>(height);

    imageWidth = width;
    imageHeight = height;
}

// ConstructRayThroughPixel is a function that returns a ray that goes through the pixel (i, j) 
Ray Config::ConstructRayThroughPixel(int i, int j, int rayNumber) {
    // Initialize variables for the top left point of the pixel and the direction of the ray.
    vec3 topLeftPoint, hitVector, rayDirection;

    // Calculate the initial top left point based on whether we're dealing with the whole pixel or sub-pixel.
    float offsetX = (rayNumber == 0) ? (pixelWidth / 2.0f) : (pixelWidth / 4.0f);
    float offsetY = (rayNumber == 0) ? (pixelHeight / 2.0f) : (pixelHeight / 4.0f);
    topLeftPoint = vec3(-1 + offsetX, 1 - offsetY, 0);

    // Calculate the hit vector based on the ray number.
    float additionalX = (rayNumber == 2 || rayNumber == 4) ? (pixelWidth / 2.0f) : 0;
    float additionalY = (rayNumber == 3 || rayNumber == 4) ? (pixelHeight / 2.0f) : 0;
    hitVector = topLeftPoint + vec3(i * pixelWidth + additionalX, -1 * (j * pixelHeight + additionalY), 0);

    // Calculate the ray direction and create the ray.
    rayDirection = normalizedVector(hitVector - eye);
    Ray ray = Ray(rayDirection, eye);

    return ray;
}

// FindIntersection is a function that returns the first intersection of the ray with the scene
Hit Config::FindIntersection(Ray ray, int ignoreObjectIndex) {
    float minT = INFINITY;
    Object* minPrimitive = new Plane(vec4(1.0, 1.0, 1.0, 1.0), Space);
    minPrimitive->objectIndex = -1;
    minPrimitive->setObjectColor(vec4(0.0, 0.0, 0.0, 0.0));
    bool gotHit = false;

    for (int i = 0; i < objects.size(); i++) {
        if (i != ignoreObjectIndex) {
            float t = objects[i]->FindIntersection(ray);
            if ((t >= 0) && (t < minT)) {
                gotHit = true;
                minPrimitive = objects[i];
                minT = t;
            }
        }
    }

    return gotHit ? Hit(ray.position + ray.direction * minT, minPrimitive)
        : Hit(ray.position + ray.direction, minPrimitive);// if there is no intersection return a hit with the space object that has the color (0,0,0,0)

}

vec4 Config::GetColor(Ray ray, Hit hit, int depth) {
    vec3 phongColor = vec3(0, 0, 0);

    // Determine the color based on the object's type
    switch (hit.object->objectType) {
    case Regular:
        phongColor = HandleRegularObject(hit, ray);
        break;
    case Reflective:
        phongColor = HandleReflectiveObject(hit, ray, depth);
        break;
    case Transparent:
        phongColor = HandleTransparentObject(hit, ray, depth);
        break;
    default:
        break;
    }

    // Return the calculated color with an alpha of 0
    return vec4(phongColor, 0.0);
}

// Handles lighting calculations for regular objects
vec3 Config::HandleRegularObject(const Hit& hit, const Ray& ray) {
    vec3 color = hit.object->getObjectColor(hit.hitPoint);
    vec3 phongModelColor = color * vec3(ambient.r, ambient.g, ambient.b); // Ambient component

    // Calculate diffuse and specular components for each light source
    for (int i = 0; i < lights.size(); i++) {
        vec3 diffuseColor = glm::max(calcDiffuseColor(hit, lights[i]), vec3(0, 0, 0));//
        vec3 specularColor = glm::max(calcSpecularColor(ray, hit, lights[i]), vec3(0, 0, 0));
        float shadow = calcShadow(hit, lights[i]);
        phongModelColor += (diffuseColor + specularColor) * shadow;
    }

    // maximum of 1.0 in each channel
    return glm::min(phongModelColor, vec3(1.0, 1.0, 1.0));
}

// Handles reflection for reflective objects
vec3 Config::HandleReflectiveObject(const Hit& hit, const Ray& ray, int depth) {
    if (depth >= 5) return vec3(0.0, 0.0, 0.0); // Recursion limit is 5

    vec3 reflectionDirection = CalculateReflectionDirection(ray.direction, hit.object->getNormal(hit.hitPoint));
    Ray reflectionRay(reflectionDirection, hit.hitPoint);
    Hit reflectedHit = FindIntersection(reflectionRay, hit.object->objectIndex);

    if (reflectedHit.object->objectType != Space) { // If the reflected ray hits an object
        vec4 reflectionColor = GetColor(reflectionRay, reflectedHit, depth + 1);
        return vec3(reflectionColor);
    }

    return vec3(0.0, 0.0, 0.0); // Default color for space
}


vec3 CalculateRefractionDirection(const vec3& incident, const vec3& normal, float eta) { 
    float cosi = dot(-incident, normal); // cosi = -N*I
    float cost2 = 1.0f - eta * eta * (1.0f - cosi * cosi); // snell's law
    if (cost2 < 0.0f) return vec3(0.0, 0.0, 0.0);// internal reflection
    return eta * incident + (eta * cosi - sqrt(cost2)) * normal; // refraction direction calculation 
}

// Checks if the ray is entering or outing the object
bool IsEntering(const vec3& incident, const vec3& normal) {
    return dot(incident, normal) < 0;
}

vec3 Config::HandleTransparentObject(const Hit& hit, const Ray& ray, int depth) {
    if (depth >= 5) return vec3(0.0, 0.0, 0.0); // Recursion limit

    float air = 1.0; 
    float material = 1.5; 
    float eta = air / material; // snell's law fraction
    vec3 normal = hit.object->getNormal(hit.hitPoint);

    // Adjust eta and normal based on ray direction
    if (!IsEntering(ray.direction, normal)) {
        std::swap(air, material);
        eta = air / material;
        normal = -normal; // Invert normal
    }

    vec3 refractedDirection = CalculateRefractionDirection(ray.direction, normal, eta);
    if (refractedDirection.length() == 0.0f) { // Handle total internal reflection
        vec3 reflectionDirection = CalculateReflectionDirection(ray.direction, normal);
        Ray reflectionRay(reflectionDirection, hit.hitPoint + reflectionDirection * 0.001f); // offset to prevent self-intersection
        Hit reflectedHit = FindIntersection(reflectionRay, hit.object->objectIndex);
        if (reflectedHit.object->objectType != Space) { // If the reflected ray hits an object
            return vec3(GetColor(reflectionRay, reflectedHit, depth + 1));
        }
    }
    else {
        // Continue with refraction
        Ray refractedRay(refractedDirection, hit.hitPoint + refractedDirection * 0.001f); // Offset to prevent self-intersection
        Hit refractedHit = FindIntersection(refractedRay, -1);
        if (refractedHit.object->objectType != Space) {
            return vec3(GetColor(refractedRay, refractedHit, depth + 1));
        }
    }
    return vec3(0.0, 0.0, 0.0); // Default color if no interaction
}


// Utility function to calculate reflection direction
vec3 Config::CalculateReflectionDirection(const vec3& incident, const vec3& normal) {
    return incident - 2.0f * normal * dot(incident, normal);
}

vec3 Config::calcDiffuseColor(Hit hit, Light* light) {
    float objectFactor = (hit.object->objectDetails.w < 0.0) ? -1.0f : 1.0f;// object factor is -1 for plane and 1 for sphere
    vec3 normalizedRayDirection = objectFactor * normalizedVector(light->lightDirection);// here we calculate the normalized ray direction based on the object factor


    // Spotlight special case
    if (light->lightType == Spot) {
        vec3 virtualSpotlightRay = normalizedVector(hit.hitPoint - light->lightPosition);
        float lightCosValue = dot(virtualSpotlightRay, objectFactor * normalizedRayDirection);

        // Checks if the spotlight hit the object
        if (lightCosValue < light->lightCosAngle) {
            return vec3(0.f, 0.f, 0.f);
        }
        normalizedRayDirection = objectFactor * virtualSpotlightRay;

    }
    vec3 objectNormal = hit.object->getNormal(hit.hitPoint);

    // N*L 
    float hitCosValue = dot(objectNormal, -normalizedRayDirection);

    // Kd*(N*L)*I_l
    vec3 diffuse_color = hit.object->getObjectColor(hit.hitPoint) * hitCosValue * light->lightRgbIntensity;
    return diffuse_color;
}

vec3 Config::calcSpecularColor(Ray ray, Hit hit, Light* light) {
    vec3 normalized_ray_direction = normalizedVector(light->lightDirection);// normalized ray direction based on the light direction

    // Spotlight
    if (light->lightType == Spot) {
        vec3 spotlight_ray = normalizedVector(hit.hitPoint - light->lightPosition); // spotlight ray is the normalized vector from the light position to the hit point
        float light_cos_value = dot(spotlight_ray, normalized_ray_direction); // L*N mulitply the normalized ray direction with spotlight ray

        // Checks if the spotlight hit the object
        if (light_cos_value < light->lightCosAngle) {// 
            return vec3(0.f, 0.f, 0.f); // Return black color - no light
        }
        normalized_ray_direction = spotlight_ray;// in case the spotlight hit the object we update the normalized ray direction to be the spotlight ray
    }
    vec3 object_normal = hit.object->getNormal(hit.hitPoint);
    vec3 reflected_light_ray =
        normalized_ray_direction - 2.0f * object_normal * dot(normalized_ray_direction, object_normal);// R
    vec3 ray_to_viewer = normalizedVector(ray.position - hit.hitPoint);//V


    float hitCosValue = dot(ray_to_viewer, reflected_light_ray);// V*R
    hitCosValue = glm::max(0.0f, hitCosValue);//alpha


    hitCosValue = pow(hitCosValue, hit.object->objectShiness);// (V*R)^n how much the object is shiny or shinniness


    float Ks = 0.7f;
    vec3 specularColor = Ks * hitCosValue * light->lightRgbIntensity;// Ks*(V*R)^n*I_l
    return specularColor;
}

float Config::calcShadow(Hit hit, Light* light) {
    vec3 normalizedRayDirection = normalizedVector(light->lightDirection);
    float minT = INFINITY;

    // Spotlight
    if (light->lightType == Spot) {
        vec3 virtualSpotlightRay = normalizedVector(hit.hitPoint - light->lightPosition);
        float lightCosValue = dot(virtualSpotlightRay, normalizedRayDirection);

        // Checks if hit the object
        if (lightCosValue < light->lightCosAngle) {
            return 0.0;
        }
        normalizedRayDirection = virtualSpotlightRay;
        // Update minT
        minT = -(dot(hit.hitPoint, light->lightPosition)) / abs(dot(-normalizedRayDirection, light->lightPosition));
    }

    // Checks if the object is in shadow
    for (int i = 0; i < objects.size(); i++) {
        if (i != hit.object->objectIndex) {
            Ray ray = Ray(-normalizedRayDirection, hit.hitPoint);
            float t = objects[i]->FindIntersection(ray);

            if ((t > 0) && (t < minT)) {
                return 0.0;
            }
        }
    }
    return 1.0;
}