// 定义了材质包：constant(单色), image(纹理映射)
#ifndef TEXTURE
#define TEXTURE

#include "vec3.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define _USE_MATH_DEFINES
#include <math.h>

class Texture
{
public:
	virtual vec3 value(double u, double v) const = 0;
};


class constant : public Texture
{
public:
	constant(const vec3& color) :color(color) {}
	vec3 value(double u, double v) const { return color; }
private:
	vec3 color;
};

class chessboard : public Texture
{
public:
	chessboard(const vec3& color, int grid) : color(color), grid(grid) {}
	vec3 value(double u, double v) const {
		u *= grid * M_PI;
		v *= grid * M_PI;
		if (sin(u) * sin(v) > 0) return color;
		else return vec3(1, 1, 1);
	}
private:
	vec3 color;
	int grid;
};


class image : public Texture
{
public:
	image(const char* fname) { img = stbi_load(fname, &w, &h, &n, 3); }
	vec3 value(double u, double v) const {
		int x = int(u * w);
		int y = int(v * h);
		if (x >= w) x = w - 1;
		if (y >= h) y = h - 1;
		return vec3(img[y*w*n + x * n],
			img[y*w*n + x * n + 1],
			img[y*w*n + x * n + 2]) / 255;
	}
	~image() { stbi_image_free(img); }
private:
	unsigned char* img; // 优先级：(h, w, n) 
	int h, w, n;
};

#endif // !TEXTURE
