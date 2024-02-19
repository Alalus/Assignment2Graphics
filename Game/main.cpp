#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "stb_image.h"
#include <iostream>
#include "../res/includes/glm/glm.hpp"
#include "Texture.h"


int main(int argc, char* argv[])
{
	const int DISPLAY_WIDTH = 512;
	const int DISPLAY_HEIGHT = 512;
	const float CAMERA_ANGLE = 0.0f;
	const float NEAR = 1.0f;
	const float FAR = 100.0f;

	Game* scn = new Game(CAMERA_ANGLE, (float)DISPLAY_WIDTH / DISPLAY_HEIGHT, NEAR, FAR);

	Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");

	Init(display);

	scn->Init();

	display.SetScene(scn);

	int width, height, numComponents;
	const char* filePath = "../res/textures/lena256.jpg";
	unsigned char* data = stbi_load(filePath, &width, &height, &numComponents, 4);
	if (data == NULL)
		std::cerr << "Unable to load texture: " << filePath << std::endl;
	unsigned char* floydFilter = floyd(data, width, height);
	unsigned char* halftoneFilter = halftone(data, width, height);
	unsigned char* cannyFilter = canny(data, width, height);

	stbi_image_free(data);

	//bottom left
	scn->AddTexture(2 * width, 2 * height, halftoneFilter);
	scn->SetShapeTex(0, 0);
	scn->CustomDraw(1, 0, scn->BACK, true, false, 0);

	//top left
	scn->AddTexture("../res/textures/lena256.jpg", false);
	scn->SetShapeTex(0, 1);
	scn->CustomDraw(1, 0, scn->BACK, false, false, 1);

	// bottom right
	scn->AddTexture(width, height, floydFilter);
	scn->SetShapeTex(0, 2);
	scn->CustomDraw(1, 0, scn->BACK, false, false, 2);

	// top right
	scn->AddTexture(width, height, cannyFilter);
	scn->SetShapeTex(0, 3);
	scn->CustomDraw(1, 0, scn->BACK, false, false, 3);

	scn->Motion();
	display.SwapBuffers();
	while (!display.CloseWindow())
	{
		//scn->Draw(1,0,scn->BACK,true,false);
		//scn->Motion();
		//display.SwapBuffers();
		display.PollEvents();

	}
	delete scn;
	// Free the memory allocated by filter functions
	delete[] floydFilter;
	delete[] halftoneFilter;
	delete[] cannyFilter;

	return 0;
}