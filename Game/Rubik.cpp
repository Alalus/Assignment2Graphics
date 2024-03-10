#include "Rubik.h"
#include <iostream>
#include "scene.h"


Rubik::Rubik()
{
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				cube[i][j][k] = 0;
			}
		}
	}
	Clockwize = 1;
	angle = 90;
}


void Rubik::setCube(int x, int y, int z, int index)
{
	cube[x][y][z] = index;

}

void Rubik::printCube()
{
	int original = 0;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				std::cout << original << "- " << "cube index: " << cube[i][j][k] << ", [" << i << "," << j << "," << k << "]" << std::endl;
				original++;
			}
		}
	}
}

void Rubik::setClockwize()
{
	Clockwize = (Clockwize == 1) ? 0 : 1;
}

void Rubik::setAngle(float angle)
{
	if (angle == 45 || angle == 90 || angle == 180)
		this->angle = angle;
	else
		std::cout << "Invalid angle" << std::endl;
}

float Rubik::getAngle()
{
	return angle;
}
glm::mat3 Rubik::rotateMatrixClockwise(glm::mat3 matrix) {


	// Transpose matrix
	for (int i = 0; i < 3; ++i) {
		for (int j = i; j < 3; ++j) {
			int temp = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = temp;
		}
	}

	// Reverse each row
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 1.5; ++j) { // Only iterate half of the row to avoid swapping twice
			int temp = matrix[i][j];
			matrix[i][j] = matrix[i][2 - j];
			matrix[i][2 - j] = temp;
		}
	}
	return matrix;
}
glm::mat3 Rubik::rotateMatrixCounterClockwise(glm::mat3 matrix) {

	// Reverse each row
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 1.5; ++j) { // Only iterate half of the row to avoid swapping twice
			int temp = matrix[i][j];
			matrix[i][j] = matrix[i][2 - j];
			matrix[i][2 - j] = temp;
		}
	}

	// Transpose matrix
	for (int i = 0; i < 3; ++i) {
		for (int j = i; j < 3; ++j) {
			int temp = matrix[i][j];
			matrix[i][j] = matrix[j][i];
			matrix[j][i] = temp;
		}
	}

	return matrix;
}

glm::mat3 Rubik::populateTempWall(int indexLayer) {
	glm::mat3 tempWall;
	switch (indexLayer) {
	case FRONT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[i][2 - j][2]; // Adjusted indexing to match the case 0 logic
		break;
	case BACK:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[2 - i][2 - j][0]; // Adjusted indexing to match the case 1 logic
		break;

	case LEFT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[0][2 - j][i]; // Adjusted indexing to match the case 1 logic
		break;
	case RIGHT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[2][2 - j][2 - i]; // Adjusted indexing to match the case 3 logic
		break;
	case UP:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[i][2][j]; // Adjusted indexing to match the case 3 logic
		break;
	case DOWN:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tempWall[j][i] = cube[i][0][2 - j]; // Adjusted indexing to match the case 3 logic
		break;
	}
	return tempWall;
}

void Rubik::populateCube(int indexLayer, glm::mat3 tempWall)
{
	switch (indexLayer)
	{
	case FRONT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				cube[i][2 - j][2] = tempWall[j][i]; // Adjusted indexing to match the case 0 logic
		break;
	case BACK:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				cube[2 - i][2 - j][0] = tempWall[j][i]; // Adjusted indexing to match the case 1 logic
		break;

	case LEFT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				cube[0][2 - j][i] = tempWall[j][i]; // Adjusted indexing to match the case 1 logic
		break;
	case RIGHT:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				cube[2][2 - j][2 - i] = tempWall[j][i]; // Adjusted indexing to match the case 3 logic
		break;
	case UP:
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				cube[i][2][j] = tempWall[j][i]; // Adjusted indexing to match the case 3 logic
			}
		}
		break;
	case DOWN:
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				cube[i][0][2 - j] = tempWall[j][i]; // Adjusted indexing to match the case 3 logic
		break;
	}
}

void Rubik::updateCube(int indexLayer)
{
	if (Clockwize)
		lockSystem.setWallAngel(lockSystem.getWallAngel(indexLayer) + angle, indexLayer);
	else
		lockSystem.setWallAngel(lockSystem.getWallAngel(indexLayer) - angle, indexLayer);
	// this is to make sure the angle is between 0 and 360
	lockSystem.modulusAngel(indexLayer);
	std::cout << "wall angle after mod: " << lockSystem.getWallAngel(indexLayer) << std::endl;
	// 
	if (lockSystem.getWallAngel(indexLayer) % 90 != 0) {
		lockSystem.lockWalls(indexLayer);
		return;
	}

	//counting the number of 90 degrees rotates
	int rotates = lockSystem.getWallAngel(indexLayer) / 90;
	rotates < 0 ? rotates = rotates * (-1) : rotates = rotates;

	glm::mat3 tempWall = populateTempWall(indexLayer);

	switch (indexLayer) {
		case FRONT:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}
			populateCube(indexLayer, tempWall);
			break;
		case BACK:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}

			populateCube(indexLayer, tempWall);
			break;
		case LEFT:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}

			populateCube(indexLayer, tempWall);
			break;

		case RIGHT:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}

			populateCube(indexLayer, tempWall);
			break;

		case UP:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}

			populateCube(indexLayer, tempWall);
			break;
		case DOWN:
			for (int i = 0; i < rotates; i++) {
				if (lockSystem.getWallAngel(indexLayer) > 0) {
					tempWall = rotateMatrixClockwise(tempWall);
				}
				else {
					tempWall = rotateMatrixCounterClockwise(tempWall);
				}
			}

			populateCube(indexLayer, tempWall);
			break;
		default:
			break;
	}
	lockSystem.unlockWalls(indexLayer);
	lockSystem.setWallAngel(0, indexLayer);
	printCube();
}

Rubik::~Rubik(void)
{
}

void Rubik::rotateLayer(int layerIndex, std::vector<Shape*> shapes) {


	float angle = getAngle();
	glm::vec3 axis;
	int index;

	// Determine rotation axis and layer based on the wall
	switch (layerIndex) {
	case FRONT:
	case BACK:
		axis = glm::vec3(0, 0, 1);
		index = (layerIndex == 0) ? 2 : 0;
		break;
	case LEFT:
	case RIGHT:
		axis = glm::vec3(1, 0, 0);
		index = (layerIndex == 2) ? 0 : 2;
		break;
	case UP:
	case DOWN:
		axis = glm::vec3(0, 1, 0);
		index = (layerIndex == 4) ? 2 : 0;
		break;
	default:
		return;
	}

	// Adjust rotation direction based on the wall
	bool isClockwise = (layerIndex == 1 || layerIndex == 2 || layerIndex == 5) ? Clockwize : !Clockwize;
	float finalAngle = isClockwise ? angle : -angle;

	// Perform rotation
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			Shape* s = nullptr;
			switch (layerIndex) {
			case 0: case 1: s = shapes[cube[i][j][index]]; break;
			case 2: case 3: s = shapes[cube[index][i][j]]; break;
			case 4: case 5: s = shapes[cube[i][index][j]]; break;
			default: return;
			}
			if (s != nullptr) {
				glm::mat4 rot_t = glm::transpose(s->getRotate());
				isClockwise ? s->MyRotate(angle, glm::vec3(rot_t * glm::vec4(axis, 1)), 0) : s->MyRotate(-angle, glm::vec3(rot_t * glm::vec4(axis, 1)), 0);
			}
		}
	}
}


void Rubik::handleR(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(RIGHT)) {
		std::cout << "cant rotate, right wall is locked" << std::endl;
		return;
	}
	rotateLayer(RIGHT, shapes);
	updateCube(RIGHT);
}
void Rubik::handleL(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(LEFT)) {
		std::cout << "cant rotate, left wall is locked" << std::endl;
		return;
	}
	rotateLayer(LEFT, shapes);
	updateCube(LEFT);

}
void Rubik::handleF(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(FRONT)) {
		std::cout << "cant rotate, front wall is locked" << std::endl;
		return;
	}
	rotateLayer(FRONT, shapes);
	updateCube(FRONT);

}
void Rubik::handleB(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(BACK)) {
		std::cout << "cant rotate, back wall is locked" << std::endl;
		return;
	}
	rotateLayer(BACK, shapes);
	updateCube(BACK);
}
void Rubik::handleU(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(UP)) {
		std::cout << "cant rotate, up wall is locked" << std::endl;
		return;
	}
	rotateLayer(UP, shapes);
	updateCube(UP);
}
void Rubik::handleD(std::vector<Shape*> shapes) {
	if (lockSystem.isWallLocked(DOWN)) {
		std::cout << "cant rotate, down wall is locked" << std::endl;
		return;
	}
	rotateLayer(DOWN, shapes);
	updateCube(DOWN);

}
void Rubik::handleSPACE() {
	setClockwize();
}
void Rubik::handleZ() {
	setAngle(angle / 2);

}
void Rubik::handleA() {
	setAngle(angle * 2);

}

void Rubik::handleRandomSequence(std::vector<Shape*> shapes) {
	for(int i = 0; i < 20; i++) {
		int random = rand() % 6; 
		switch (random) {
		case 0: handleR(shapes); break;
		case 1: handleL(shapes); break;
		case 2: handleF(shapes); break;
		case 3: handleB(shapes); break;
		case 4: handleU(shapes); break;
		case 5: handleD(shapes); break;
		}
	}
}
