#include "InputManager.h"
#include "game.h"
#include "glm/glm.hpp"
#include "Config.h"
#include "Ray.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;
using namespace glm;
int main(int argc, char* argv[])
{
    const int DISPLAY_WIDTH = 800;
    const int DISPLAY_HEIGHT = 800;
    const float CAMERA_ANGLE = 0.0f;
    const float NEAR = 1.0f;
    const float FAR = 100.0f;

    Game* scn = new Game(CAMERA_ANGLE, (float)DISPLAY_WIDTH / DISPLAY_HEIGHT, NEAR, FAR);

    Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");

    Init(display);

    scn->Init();

    display.SetScene(scn);
    const char* fileName = "../res/txtfiles/custom_scene.txt";
    //------------------------- step 1: read the scene file and create the scene -------------------------//

    Config sceneConfigure = Config();
    // read the scene file and create the scene
    sceneConfigure.readSceneFile(fileName, DISPLAY_WIDTH, DISPLAY_HEIGHT);

    //--------------------------------------------------------------------------------------------------//


    //------------------------- step 2: render the scene -------------------------//

    int imageWidth = sceneConfigure.imageWidth, imageHeight = sceneConfigure.imageHeight;
    // create a matrix to store the color of each pixel
    std::vector<std::vector<vec4>> pixelColorMatrix(imageHeight, std::vector<vec4>(imageWidth));

    // run through each pixel and calculate the color of the pixel 
    for (int j = 0; j < imageHeight; j++) {
        for (int i = 0; i < imageWidth; i++) {

            vec4 pixelColor;

            if (sceneConfigure.bonusFlag < 1.0) {
                // ConstructRayThroughPixel is a function that returns a ray that goes through the pixel (i, j)
                Ray ray = sceneConfigure.ConstructRayThroughPixel(i, j, 0);

                // FindIntersection is a function that returns the first intersection of the ray with the scene
                Hit hit = sceneConfigure.FindIntersection(ray, -1);

                // GetColor is a function that returns the color of the pixel (i, j)
                pixelColorMatrix[j][i] = sceneConfigure.GetColor(ray, hit, 0);
            }
            //Bonus
            else {
                for (int currentExtraRay = 1; currentExtraRay < 5; currentExtraRay++) {
                    Ray currentRay = sceneConfigure.ConstructRayThroughPixel(i, j, currentExtraRay);
                    Hit currentHit = sceneConfigure.FindIntersection(currentRay, -1);
                    vec4 currentPixelColor = sceneConfigure.GetColor(currentRay, currentHit, 0);
                    pixelColor += currentPixelColor;
                }
                pixelColor = pixelColor / 4.0f; // average the color of the 4 rays
            }
        }
    }
    // create an image from the pixelColorMatrix and add it to the scene as a texture 



    unsigned char* image = new unsigned char[imageWidth * imageHeight * 4];
    int index = 0;
    // fill the image array with the pixelColorMatrix values
    for (int i = 0; i < imageHeight; i++) {
        for (int j = 0; j < imageWidth; j++) {
            image[index] = (unsigned char)(pixelColorMatrix[i][j].x * 255);
            image[index + 1] = (unsigned char)(pixelColorMatrix[i][j].y * 255);
            image[index + 2] = (unsigned char)(pixelColorMatrix[i][j].z * 255);
            image[index + 3] = (unsigned char)(pixelColorMatrix[i][j].w * 255);
            index += 4;
        }

    }

    // add the image to the scene as a texture
    scn->AddTexture(imageWidth, imageHeight, image);

    while (!display.CloseWindow())
    {
        scn->Draw(1, 0, scn->BACK, true, false);
        scn->Motion();
        display.SwapBuffers();
        display.PollEvents();

    }
    delete scn;
    return 0;
}