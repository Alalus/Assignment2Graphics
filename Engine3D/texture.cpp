
#define STB_IMAGE_IMPLEMENTATION
#include "texture.h"
#include "stb_image.h"
#include "../res/includes/glad/include/glad/glad.h"
#include <iostream>
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <numbers>
#include <fstream>

void rgbaToGrey(std::vector<std::vector<unsigned char>>& greyScaleImg, unsigned char** rgbImg, int square_size);
void greyToRgba(unsigned char** rgbImg, std::vector<std::vector<unsigned char>>& greyScaleImg, unsigned char** data, int square_size);
int findDirection(float dir);
int findNearestColor(int originalIntensity);

// Assuming the input data is RGBA (4 channels), convert it to grayscale
void rgbaToGrey(std::vector<std::vector<unsigned char>>& greyScaleImg, unsigned char** rgbImg, int square_size) {
    for (int y = 0; y < square_size; y++) {
        for (int x = 0; x < square_size; x++) {
            // Calculate the correct index for the 1D RGBA data array
            int index = (y * square_size + x) * 4; // 4 channels per pixel
            // Compute the grayscale value as the average of the RGB components
            unsigned char grayscaleValue = static_cast<unsigned char>
                (((*rgbImg)[index] + (*rgbImg)[index + 1] + (*rgbImg)[index + 2]) / 3);
            // Assign the grayscale value to the matrix
            greyScaleImg[y][x] = grayscaleValue;
        }
    }
}
//Convert grayscale back to RGBA
void greyToRgba(unsigned char** rgbImg, std::vector<std::vector<unsigned char>>& greyScaleImg, unsigned char** data, int square_size) {
    int index = 0;
    for (int i = 0; i < square_size; i++) {
        for (int j = 0; j < square_size; j++) {
            //red value
            (*rgbImg)[index] = greyScaleImg[i][j];
            //green value
            (*rgbImg)[index + 1] = greyScaleImg[i][j];
            //blue value
            (*rgbImg)[index + 2] = greyScaleImg[i][j];
            //alpha value
            (*rgbImg)[index + 3] = ((*data)[4 * (((i / 2) * (square_size / 2)) + (j / 2)) + 3]);
            index += 4;
        }
    }
}
// Find the direction of the gradient for Canny algorithm 
int findDirection(float dir) {
    int caseNumber;
    dir *= (180.0 / M_PI); // Convert from radians to degrees
    float teta = 180.0 / 8;
    // Normalize the direction to be within [0, 360) degrees
    dir = fmod(dir + 360, 360.0);
    if ((0 <= dir && dir < teta) || (315 + teta <= dir && dir <= 360) || (135 + teta <= dir && dir < 180 + teta)) { // horizontal edge
        caseNumber = 0;
    }
    else if ((teta <= dir && dir < 90 - teta) || (180 + teta <= dir && dir < 225 + teta)) { // 45 degree right diagonal
        caseNumber = 1;
    }
    else if ((45 + teta <= dir && dir < 90 + teta) || (225 + teta <= dir && dir < 270 + teta)) { // vertical edge
        caseNumber = 2;
    }
    else { // 45 degree left diagonal
        caseNumber = 3; 
    }
    return caseNumber;
}
// Find the nearest color for Floyd-Steinberg Dithering
int findNearestColor(int originalIntensity) {
    return static_cast<int>((originalIntensity / 16) * 16);
}

//=======================================================================================================ALGORITHMS=======================================================================================================

// Perform floyd on the input image data
unsigned char* floyd(unsigned char* data, int width, int height) {
    // Initialize a 2D vector for the grayscale values
    std::vector<std::vector<unsigned char>> grayscaleMatrix(height, std::vector<unsigned char>(width));
    rgbaToGrey(grayscaleMatrix, &data, width);

    // Apply Floyd-Steinberg Dithering
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            unsigned char originalPixel = grayscaleMatrix[i][j];
            unsigned char newPixel = findNearestColor(originalPixel);
            grayscaleMatrix[i][j] = newPixel;
            int err = originalPixel - newPixel;

            if (j + 1 < width) // not the last element in row
                grayscaleMatrix[i][j + 1] = std::min(255, std::max(0, grayscaleMatrix[i][j + 1] + err * 7 / 16));
            if (i + 1 < height) { // not the last row
                if (j > 0)
                    grayscaleMatrix[i + 1][j - 1] = std::min(255, std::max(0, grayscaleMatrix[i + 1][j - 1] + err * 3 / 16)); 
                grayscaleMatrix[i + 1][j] = std::min(255, std::max(0, grayscaleMatrix[i + 1][j] + err * 5 / 16)); 
                if (j + 1 < width)
                    grayscaleMatrix[i + 1][j + 1] = std::min(255, std::max(0, grayscaleMatrix[i + 1][j + 1] + err * 1 / 16));
            }
        }
    }

    // Write the pixel values to a file
    std::ofstream outputFile("img6.txt");
    if (outputFile.is_open()) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                outputFile << static_cast<int>(grayscaleMatrix[i][j] / 16);
                if (j < width - 1) outputFile << ",";
            }
            outputFile << "\n";
        }
        outputFile.close();
    }
    else {
        std::cerr << "Unable to open file for writing." << std::endl;
        return nullptr;
    }

    // Convert the grayscale matrix back to RGBA
    unsigned char* outputData = new unsigned char[width * height * 4];
    greyToRgba(&outputData, grayscaleMatrix, &data, width);

    return outputData; // Remember to manage the memory for outputData
}
// Perform halftoning on the input image data
unsigned char* halftone(unsigned char* data, int width, int height) {

    //Create a matrix to store grayscale values of the original image
    std::vector<std::vector<unsigned char>> grayscaleMatrix(height, std::vector<unsigned char>(width, 0));
    rgbaToGrey(grayscaleMatrix, &data, width);

    // Create a matrix to store grayscale values after halftoning
    std::vector<std::vector<unsigned char>> imgAfterHalftone(2 * height, std::vector<unsigned char>(2 * width, 0));

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float normalizedPixelValue = (float)grayscaleMatrix[i][j] / (float)255;
            i *= 2;
            j *= 2;
            // Halftoning logic based on pixel intensity
            if (normalizedPixelValue < 0.2) {
                imgAfterHalftone[i][j] = 0;
                imgAfterHalftone[i + 1][j] = 0;
                imgAfterHalftone[i][j+ 1] = 0;
                imgAfterHalftone[i + 1][j + 1] = 0;
            }
            else if (normalizedPixelValue < 0.4) {
                imgAfterHalftone[i][j] = 0;
                imgAfterHalftone[i + 1][j] = 255;
                imgAfterHalftone[i][j + 1] = 0;
                imgAfterHalftone[i + 1][j + 1] = 0;

            }
            else if (normalizedPixelValue < 0.6) {
                imgAfterHalftone[i][j] = 0;
                imgAfterHalftone[i + 1][j] = 255;
                imgAfterHalftone[i][j + 1] = 255;
                imgAfterHalftone[i + 1][j + 1] = 0;

            }
            else if (normalizedPixelValue < 0.8) {
                imgAfterHalftone[i][j] = 0;
                imgAfterHalftone[i + 1][j] = 255;
                imgAfterHalftone[i][j + 1] = 255;
                imgAfterHalftone[i + 1][j + 1] = 255;

            }
            else {
                imgAfterHalftone[i][j] = 255;
                imgAfterHalftone[i + 1][j] = 255;
                imgAfterHalftone[i][j + 1] = 255;
                imgAfterHalftone[i + 1][j + 1] = 255;
            }
            i /= 2;
            j /= 2;
        }
    }
    // Write the pixel values to a file
    std::ofstream outputFile("img5.txt");
    if (outputFile.is_open()) {
        for (int i = 0; i < height*2; i++) {
            for (int j = 0; j < width*2; j++) {
                outputFile << static_cast<int>(imgAfterHalftone[i][j] / 255);
                if (j < width*2 - 1) outputFile << ",";
            }
            outputFile << "\n";
        }
        outputFile.close();
    }
    else {
        std::cerr << "Unable to open file for writing." << std::endl;
        return nullptr;
    }
    unsigned char* halftoneOutput = new unsigned char[width * height * 16];

    // Convert the halftoned grayscale image back to RGBA
    greyToRgba(&halftoneOutput, imgAfterHalftone, &data, width * 2);

    return halftoneOutput;
}
// Perform canny on the input image data
unsigned char* canny(unsigned char* data, int width, int height) {

    std::vector<std::vector<unsigned char>> grayscaleMatrix(height, std::vector<unsigned char>(width, 0));
    rgbaToGrey(grayscaleMatrix, &data, width); // Convert the input image to grayscale

    std::vector<std::vector<unsigned char>> blurredMatrix(height, std::vector<unsigned char>(width, 0));

    // Gaussian filter
    float gaussKernel[3][3] = {
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };
    float kernelSum = 16.0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float filteredPixel = 0.0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int posY = y + dy;
                    int posX = x + dx;
                    // Handle the borders
                    posY = std::max(0, std::min(height - 1, posY));
                    posX = std::max(0, std::min(width - 1, posX));
                    filteredPixel += grayscaleMatrix[posY][posX] * gaussKernel[dy + 1][dx + 1]; 
                }
            }
            // Normalize and assign the new pixel value 
            blurredMatrix[y][x] = static_cast<unsigned char>(filteredPixel / kernelSum);
        }
    }

    float Gx[3][3] = { {-1, 0, 1},
                      {-2, 0, 2},
                      {-1, 0, 1} };
    float Gy[3][3] = { {1, 2, 1},
                      {0, 0, 0},
                      {-1, -2, -1} };

    std::vector<std::vector<float>> gradientMagnitude(height, std::vector<float>(width, 0));
    std::vector<std::vector<float>> gradientDirection(height, std::vector<float>(width, 0));

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float sumX = 0.0;
            float sumY = 0.0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int posY = y + dy;
                    int posX = x + dx;
                    // Handle the borders
                    posY = std::max(0, std::min(height - 1, posY));
                    posX = std::max(0, std::min(width - 1, posX));
                    sumX += Gx[dy + 1][dx + 1] * blurredMatrix[posY][posX];
                    sumY += Gy[dy + 1][dx + 1] * blurredMatrix[posY][posX];
                }
            }
            gradientMagnitude[y][x] = sqrt(sumX * sumX + sumY * sumY);
            gradientDirection[y][x] = atan2(sumY, sumX);
        }
    }

    std::vector<std::vector<unsigned char>> edgesMatrix(height, std::vector<unsigned char>(width, 0));

    // Non-maximum suppression
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            float dir = gradientDirection[y][x];
            float mag = gradientMagnitude[y][x];
            float p1, p2;
            int caseNumber = findDirection(dir);
            switch (caseNumber) {
            case 0: // horizontal edge
                p1 = gradientMagnitude[y][x + 1];
                p2 = gradientMagnitude[y][x - 1];
                break;
            case 1: // 45 degree right diagonal
                p1 = gradientMagnitude[y + 1][x - 1];
                p2 = gradientMagnitude[y - 1][x + 1];
                break;
            case 2: // vertical edge
                p1 = gradientMagnitude[y - 1][x];
                p2 = gradientMagnitude[y + 1][x];
                break;
            case 3: // 45 degree left diagonal
                p1 = gradientMagnitude[y + 1][x + 1];
                p2 = gradientMagnitude[y - 1][x - 1];
                break;

            }

            if (mag < p1 || mag < p2) {
                edgesMatrix[y][x] = 0; // supress
            }
            else {
                edgesMatrix[y][x] = static_cast<unsigned char>(mag);
            }
        }
    }

    // Double thresholding
    float highThres = 180;
    float lowThres = 120;

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            float mag = edgesMatrix[y][x];
            if (highThres <= mag) { // strong edge
                edgesMatrix[y][x] = 255;
            }
            else if (lowThres <= mag) {  // weak edge
                edgesMatrix[y][x] = 1; // 1 for indicate weak edge
            }
            else { // non edge
                edgesMatrix[y][x] = 0;
            }
        }
    }

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int value = edgesMatrix[y][x];
            if (value == 1) { // weak edge - check if one of 8 neighbors is strong edge
                if (edgesMatrix[y - 1][x - 1] == 255 || edgesMatrix[y - 1][x] == 255 || edgesMatrix[y - 1][x + 1] == 255 ||
                    edgesMatrix[y][x - 1] == 255 || edgesMatrix[y][x + 1] == 255 ||
                    edgesMatrix[y + 1][x - 1] == 255 || edgesMatrix[y + 1][x] == 255 || edgesMatrix[y + 1][x + 1] == 255) {
                    edgesMatrix[y][x] = 255;
                }
                else { // suppress
                    edgesMatrix[y][x] = 0;
                }
            }
        }
    }
    // Write the pixel values to a file
    std::ofstream outputFile("img4.txt");
    if (outputFile.is_open()) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                outputFile << static_cast<int>(edgesMatrix[i][j] / 255);
                if (j < width - 1) outputFile << ",";
            }
            outputFile << "\n";
        }
        outputFile.close();
    }
    else {
        std::cerr << "Unable to open file for writing." << std::endl;
        return nullptr;
    }

    unsigned char* canny_output = new unsigned char[width * height * 4];
    greyToRgba(&canny_output, edgesMatrix, &data, width);

    return canny_output;

}

//=======================================================================================================ALGORITHMS=======================================================================================================

Texture::Texture(const std::string& fileName)
{
    int width, height, numComponents;
    unsigned char* data = stbi_load((fileName).c_str(), &width, &height, &numComponents, 4);

    if (data == NULL)
        std::cerr << "Unable to load texture: " << fileName << std::endl;

    glGenTextures(1, &m_texture);
    Bind(m_texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_LOD_BIAS, -0.4f);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    stbi_image_free(data);
}

Texture::Texture(int width, int height, unsigned char* data)
{
    glGenTextures(1, &m_texture);
    Bind(m_texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

}

Texture::~Texture()
{
    glDeleteTextures(1, &m_texture);
}

void Texture::Bind(int slot)
{
    glActiveTexture(GL_TEXTURE0 + slot);
    glBindTexture(GL_TEXTURE_2D, m_texture);
}