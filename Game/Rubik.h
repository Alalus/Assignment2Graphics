# pragma once
#include <vector>
#include "shape.h"
#include "LockSystem.h"

// Define an enumeration for the cube faces

class Rubik {
public:
	Rubik();
	void setCube(int x, int y, int z, int index);
	void printCube();
	void setClockwize();
	void setAngle(float angle);
	float getAngle();
	glm::mat3 populateTempWall(int indexLayer);
	void populateCube(int indexLayer, glm::mat3 tempWall);
	glm::mat3 rotateMatrixClockwise(glm::mat3 matrix);
	glm::mat3 rotateMatrixCounterClockwise(glm::mat3 matrix);
	void Rubik::updateCube(int indexLayer);
	~Rubik(void);
	void rotateLayer(int layerIndex, std::vector<Shape*> shapes);
	void handleR(std::vector<Shape*> shapes);
	void handleL(std::vector<Shape*> shapes);
	void handleF(std::vector<Shape*> shapes);
	void handleB(std::vector<Shape*> shapes);
	void handleU(std::vector<Shape*> shapes);
	void handleD(std::vector<Shape*> shapes);
	void handleSPACE();
	void handleZ();
	void handleA();
	void handleRandomSequence(std::vector<Shape*> shapes);

private:
	int cube[3][3][3];
	int Clockwize;
	float angle;
	LockSystem lockSystem;

};