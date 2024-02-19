#pragma once
#ifndef TEXTURE_H
#define TEXTURE_H
#include <vector>
#include <string>
unsigned char* floyd(unsigned char* data, int width, int height);
unsigned char* halftone(unsigned char* data, int width, int height);
unsigned char* canny(unsigned char* data, int width, int height);


class Texture
{
public:
	Texture(const std::string& fileName);
	Texture(int width, int height, unsigned char* data);
	void Bind(int slot);
	inline int GetSlot() { return m_texture; }
	~Texture();
protected:
private:
	Texture(const Texture& texture) {}
	void operator=(const Texture& texture) {}
	unsigned int m_texture;
	std::vector<std::vector<unsigned char>> mat;
};

#endif