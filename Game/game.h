#pragma once
#include "scene.h"
#include "Rubik.h"

class Game : public Scene
{
public:

	Game();
	Game(float angle, float relationWH, float near, float far);
	void Init();
	void Update(const glm::mat4& MVP, const glm::mat4& Model, const int  shaderIndx);
	void ControlPointUpdate();
	void WhenRotate();
	void WhenTranslate();
	void Motion();
	~Game(void);
	Rubik* GetCube() const { return cube; }

protected:

private:
	Rubik* cube = nullptr;
};