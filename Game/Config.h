#pragma once

#include <string>
#include <vector>
#include "Light.h"
#include "glm/glm.hpp"
#include "Ray.h"
#include "Object.h"
#include "Light.h"
#include "Hit.h"

using namespace std;
using namespace glm;

class Config {
public:
    //scene file reading and processing
    void Config::readSceneFile(const std::string& fileName, int width, int height);
    void Config::processSceneData(const std::vector<std::vector<std::string>>& sceneData);
    void Config::processLineType(const std::string& lineType, const vec4& dataVector);

    //scene setup methods
    void Config::addLight(const vec4& dataVector);
    void Config::addObject(const vec4& dataVector, objectTypes type);
    void Config::setColorForObjects();
    void Config::setSpotlightPositions();
    void Config::initializeScene(int width, int height);

    // Ray tracing core functionalities
    Ray ConstructRayThroughPixel(int i, int j, int rayNumber);
    Hit FindIntersection(Ray ray, int ignoreObjectIndex);
    vec4 GetColor(Ray ray, Hit hit, int depth);
    vec3 HandleRegularObject(const Hit& hit, const Ray& ray);
    vec3 HandleReflectiveObject(const Hit& hit, const Ray& ray, int depth);
    vec3 HandleTransparentObject(const Hit& hit, const Ray& ray, int depth);
    vec3 CalculateReflectionDirection(const vec3& incident, const vec3& normal);


    // Lighting calculations methods
    vec3 calcDiffuseColor(Hit hit, Light* light);
    vec3 calcSpecularColor(Ray ray, Hit hit, Light* light);
    float calcShadow(Hit hit, Light* light);

    //scene data
    vector<vector<string>> sceneData;

    // scene elements
    vec3 eye;
    float bonusFlag;

    //lighting properties
    vec4 ambient;
    vector<Object*> objects;
    vector<vec4> colors;
    vector<Light*> lights;
    vector<vec4> positions;
    vector<vec4> intensities;

    //image rendering properties
    int imageWidth, imageHeight;
    float pixelWidth, pixelHeight;
};
