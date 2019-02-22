// Meterial：决定反射/折射光的方向（光源是一种特殊的材质：不会发出光线）
// fuzzy -> diffuse / specular -> reflection / transparent -> reflection or refraction
#ifndef MATERIAL
#define MATERIAL

#include <stdlib.h>
#include <time.h>
#include "vec3.h"
#include "ray.h"
#include "texture.h"

class Material;

// 存放各种模型属性，用做求交点后的返回参数包
struct hitParam
{
	double t;
	vec3 pos;
	vec3 normal;
	double u, v;
	Texture* tp;
	Material* mp;
};

class Material
{
public:
	virtual bool scatter(const ray& in, const hitParam& param, ray& out) const = 0;
	virtual bool is_light() const = 0;
};

// light 是一种特殊的材质=。=别的材质没有scatter出ray的时候是黑色，但light是光的颜色
class light : public Material
{
public:
	bool scatter(const ray& in, const hitParam& param, ray& out) const { return false; }
	bool is_light() const { return true; }
};

class fuzzy : public Material
{
public:
	fuzzy() { srand(unsigned(time(NULL))); }
	bool scatter(const ray& in, const hitParam& param, ray& out) const;
	bool is_light() const { return false; }
};

class specular : public Material
{
public:
	specular(double fuzz): fuzz(fuzz) { srand(unsigned(time(nullptr))); }
	bool scatter(const ray& in, const hitParam& param, ray& out) const;
	bool is_light() const { return false; }
private:
	double fuzz;
};

class transparent : public Material
{
public:
	transparent(double ref) :ref_idx(ref) {}
	bool scatter(const ray& in, const hitParam& param, ray& out) const;
	bool is_light() const { return false; }
private:
	double ref_idx;
};

/*************************** implementation **************************/
// 在半径为1的球中随机采样
vec3 randam_vector() {
	vec3 p;
	double x, y, z;
	do {
		x = rand() / double(RAND_MAX);
		y = rand() / double(RAND_MAX);
		z = rand() / double(RAND_MAX);
		p = vec3(x, y, z) * 2 - vec3(1, 1, 1);
	} while (p.length() >= 1);
	return p;
}

inline vec3 reflect(const vec3& in, const vec3& normal) {
	double n_length = -in.dot(normal) / normal.length();
	return in + normal * n_length * 2;
}

bool refract(const vec3& in, const vec3& normal, double ni_over_nr, vec3& out) {
	double cosInTheta = -normalize(in).dot(normal);
	double sinOutTheta2 = ni_over_nr * ni_over_nr * (1 - cosInTheta * cosInTheta);
	if (sinOutTheta2 < 1) {
		out = -normal * sqrt(1 - sinOutTheta2) + (normal * cosInTheta + normalize(in)) * ni_over_nr;
		return true;
	}
	else return false;
}

// 用于计算折射/反射概率
double schlick(double cos, double ref_idx) {
	double r0 = (1 - ref_idx) / (1 + ref_idx);
	r0 = r0 * r0;
	return r0 + (1 - r0)*pow((1 - cos), 5);
}

/************************* fuzzy *********************************/
bool fuzzy::scatter(const ray& in, const hitParam& param, ray& r) const {
	vec3 diffuse = param.pos + param.normal + randam_vector();
	r = ray(param.pos, diffuse - param.pos);
	return true;
}

/************************ specular *********************************/
bool specular::scatter(const ray& in, const hitParam& param, ray& out) const {
	double max_r = -in.direction().dot(param.normal);
	if (max_r <= 0) return false;
	vec3 out_direction = reflect(in.direction(), param.normal) + randam_vector() * max_r * fuzz;
	out = ray(param.pos, out_direction);
	return true;
}

/************************ transparent ******************************/
bool transparent::scatter(const ray& in, const hitParam& param, ray& out) const {
	vec3 reflection = reflect(in.direction(), param.normal);
	double delt = in.direction().dot(param.normal);
	vec3 outward_normal = param.normal;
	double ni_over_nr = 1 / ref_idx;
	double cos = -delt / in.direction().length();
	if (delt > 0) {
		outward_normal = -outward_normal;
		ni_over_nr = ref_idx;
		cos = sqrt(1 - ref_idx * ref_idx*(1 - cos * cos));
	}
	vec3 refraction;
	if (refract(in.direction(), outward_normal, ni_over_nr, refraction)) {
		double prob = schlick(cos, ref_idx);
		// double prob = 1;
		if (double(rand()) / RAND_MAX < prob) {
			out = ray(param.pos, reflection);
			return true;
		}
	}
	out = ray(param.pos, refraction);
	return true;
}

#endif