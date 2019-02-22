// 定义了各种模型，以及如何计算交点
#ifndef OBJECT
#define OBJECT

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include "vec3.h"
#include "ray.h"
#include "material.h"
#include "texture.h"

using namespace std;

class Object
{
public:
	virtual bool hit(const ray& in, double min, double max, hitParam& param) const = 0;
	virtual vec3 center() const = 0;
};


class Object_list
{
public:
	Object_list() {};
	void add(Object* op) { list.push_back(op); }
	bool hit(const ray& in, double min, double max, hitParam& param, int& idx) const {
		double current_max = max;
		bool flag = false;
		for (int i = 0; i < list.size(); ++i) {
			if (list[i]->hit(in, min, current_max, param)) {
				current_max = param.t;
				idx = i;
				flag = true;
			}
		}
		return flag;
	}
	Object*& operator[] (unsigned i) { return list[i]; }
private:
	vector<Object*> list;
};


class Sphere: public Object
{
public:
	Sphere(const vec3& origin, double radius, Material* m, Texture* t) :origin(origin), r(radius), mp(m), tp(t) {};
	vec3 center() const { return origin; }
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		double a = in.direction().squared_length();
		double b = 2 * in.direction().dot(in.start() - origin);
		double c = (in.start() - origin).squared_length() - r * r;
		double delta = b * b - 4 * a*c;
		if (delta <= 0) return false;
		double t = (-b - sqrt(delta)) / 2 / a;
		if (t < min || t > max) {
			t = (-b + sqrt(delta)) / 2 / a;
			if (t < min || t > max)
				return false;
		}
		param.t = t;
		param.pos = in.position(t);
		param.normal = normalize(param.pos - origin);
		param.mp = mp;
		param.tp = tp;
		get_uv(param.pos, param.u, param.v);
		return true;
	}
	void get_uv(vec3 p, double& u, double& v) const {
		vec3 unit_p = normalize(p - origin);
		double theta = acos(unit_p[1]);		// [0, PI]
		double alpha = atan2(-unit_p[2], unit_p[0]); // [-PI, PI]
		v = theta / M_PI;
		u = (alpha + M_PI) / (M_PI * 2);
	}
private:
	vec3 origin;
	double r;
	Material* mp;
	Texture* tp;
};


class rect : public Object
{
public:
	rect(int align, double k, double u0, double u1, double v0, double v1, Material* m, Texture* t)
		:align(align), k(k), u0(u0), u1(u1), v0(v0), v1(v1), mp(m), tp(t) {
		for(int i = 0;i < 3;++i)
			if (i != align) { i0 = i; break; }
		for (int i = 0; i < 3; ++i)
			if (i != align && i != i0) i1 = i;
		normal = vec3(0, 0, 0);
		normal[align] = 1;
	}
	vec3 center() const {
		vec3 center = vec3(0, 0, 0);
		center[align] = k;
		center[i0] = (u0 + u1) / 2;
		center[i1] = (v0 + v1) / 2;
		return center;
	}
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		// 计算和平面的交点
		double t = (k - in.start()[align]) / in.direction()[align];
		if (t <= min || t >= max)
			return false;
		vec3 p = in.position(t);
		// 判断交点是否在矩形内
		if (p[i0] >= u0 && p[i0] <= u1 && p[i1] >= v0 && p[i1] <= v1) {
			param.t = t;
			param.pos = p;
			param.mp = mp;
			param.tp = tp;
			param.normal = normal;
			param.u = (p[i0] - u0) / (u1 - u0);
			param.v = (v1 - p[i1]) / (v1 - v0);
			return true;
		}
		return false;
	}
private:
	int align;
	int i0, i1;
	double k;
	double u0, u1, v0, v1;
	vec3 normal;  // 法线默认坐标轴正方向
	Material *mp;
	Texture *tp;
};


class cube : public Object {
public:
	cube(vec3 origin, double l, Material* m, Texture *t) :origin(origin), length(l), mp(m), tp(t) {
		face.add(new rect(2, origin[2] + l, origin[0], origin[0] + l, origin[1], origin[1] + l, m, t));  // front
		face.add(new rect(2, origin[2], origin[0], origin[0] + l, origin[1], origin[1] + l, m, t));      // back
		face.add(new rect(0, origin[0], origin[1], origin[1] + l, origin[2], origin[2] + l, m, t));		 // left
		face.add(new rect(0, origin[0] + l, origin[1], origin[1] + l, origin[2], origin[2] + l, m, t));	 // right
		face.add(new rect(1, origin[1], origin[0], origin[0] + l, origin[2], origin[2] + l, m, t));		 // bottom
		face.add(new rect(1, origin[1] + l, origin[0], origin[0] + l, origin[2], origin[2] + l, m, t));	 // up
	}
	vec3 center() const { return origin + vec3(length/2, length/2, length/2); }
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		int idx;
		if (face.hit(in, min, max, param, idx)) {
			// 翻转法向量，因为有一半的面法向量应该朝负半轴方向
			if (idx == 1 || idx == 2 || idx == 4)
				param.normal = -param.normal;
			// get_uv(idx, param.u, param.v);
			return true;
		}
		return false;
	}
private:
	vec3 origin;
	double length;
	Object_list face;
	Material* mp;
	Texture* tp;

	void get_uv(int idx, double& u, double& v) const {
		// 0 1 2
		// 3 4 5
		double u_offset = (idx % 3) / 3.0;
		double v_offset = 0.5 * (idx / 3);
		u = u_offset + u / 3;
		v = v_offset + v / 2;
	}
};


class triangle : public Object
{
public:
	triangle(const vec3& p1, const vec3& p2, const vec3& p3, const vec3& vn1, const vec3& vn2, const vec3& vn3,
		double u1, double u2, double u3, double v1, double v2, double v3)
		:p1(p1), p1p2(p2-p1), p1p3(p3-p1), n1(vn1), n2(vn2), n3(vn3), u1(u1), u2(u2), u3(u3), v1(v1), v2(v2), v3(v3) {
			normal = p1p2.cross(p1p3);
	}
	vec3 center() const { return p1 + (p1p2 + p1p3) / 3; }
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		vec3 p1s = in.start() - p1;
		vec3 P = in.direction().cross(p1p3);
		vec3 Q = p1s.cross(p1p2);
		double delt = P.dot(p1p2);
		// delt < 0 ray从背面入射，|delt| < 1e-10 无解
		if (delt < 1e-10) return false;

		double invDelt = 1 / delt;
		// t: 判断是否与平面相交
		double t = Q.dot(p1p3) * invDelt;
		if (t < min || t > max) return false;
		// u: 判断是否在三角形内
		double u = P.dot(p1s) * invDelt;
		if (u < 0 || u > 1) return false;
		// v: 判断是否在三角形内
		double v = Q.dot(in.direction()) * invDelt;
		if (v < 0 || u + v > 1) return false;
		double alpha = 1 - u - v;
		param.t = t;
		param.pos = in.position(t);
		param.normal = n1 * alpha + n2 * u + n3 * v;
		param.u = alpha * u1 + u * u2 + v * u3;
		param.v = alpha * v1 + u * v2 + v * v3;

		return true;
	}
private:
	vec3 p1, p1p2, p1p3;
	vec3 n1, n2, n3;
	double u1, u2, u3;
	double v1, v2, v3;
	vec3 normal;
};


class model : public Object
{
public:
	model(const char* fname, int scaleRate, Material* m, Texture* t) : mp(m), tp(t) {
		ifstream f(fname);
		if (f.is_open()) {
			string tmp;
			vector<double> u, v;
			vector<vec3> pos, vn;
			
			// skip header
			while ((f >> tmp) && tmp != "v");
			// read vertices position
			double x, y, z;
			do {
				f >> x >> y >> z;
				pos.push_back(vec3(x, y, z) * scaleRate);
			} while ((f >> tmp) && tmp == "v");
			// read texture coordinates
			double tu, tv;
			do {
				f >> tu >> tv;
				u.push_back(tu);
				v.push_back(tv);
			} while ((f >> tmp) && tmp == "vt");
			// read normal
			do {
				f >> x >> y >> z;
				vn.push_back(vec3(x, y, z));
			} while ((f >> tmp) && tmp == "vn");
			// skip mtl
			while ((f >> tmp) && tmp != "f");
			// read meshes
			int pi[3], ti[3], ni[3];
			do {
				for (int i = 0; i < 3; ++i) {
					f >> tmp;
					sscanf_s(tmp.c_str(), "%d/%d/%d/", &pi[i], &ti[i], &ni[i]);
				}
				mesh.add(new triangle(pos[pi[0]-1], pos[pi[1]-1], pos[pi[2]-1], vn[ni[0]-1], vn[ni[1]-1], vn[ni[2]-1],
					u[ti[0]-1], u[ti[1]-1], u[ti[2]-1], v[ti[0]-1], v[ti[1]-1], v[ti[2]-1]));
			} while ((f >> tmp) && tmp == "f");
			// build bounding box
			vec3 min, max;
			get_min_max_avg(pos, min, max, avg);
			vec3 dxyz = max - min;
			bbox = new cube(min, dxyz.max(), NULL, NULL);
		}
	}
	vec3 center() const { return avg; }
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		if (bbox->hit(in, min, max, param)) {
			int _;
			if (mesh.hit(in, min, max, param, _)) {
				param.v = 1 - param.v;		// obj文件里v是向上的
				param.mp = mp;
				param.tp = tp;
				return true;
			}
		}
		return false;
	}

private:
	cube* bbox;
	Object_list mesh;
	vec3 avg;
	Material* mp;
	Texture* tp;

	void get_min_max_avg(const vector<vec3>& pos, vec3& min, vec3& max, vec3& avg) const {
		min = max = pos[0];
		avg = vec3(0, 0, 0);
		for (int i = 0; i < pos.size(); ++i) {
			avg = avg + pos[i];
			for (int j = 0; j < 3; ++j) {
				if (pos[i][j] < min[j]) min[j] = pos[i][j];
				else if (pos[i][j] > max[j]) max[j] = pos[i][j];
			}
		}
		avg = avg / pos.size();
	}
};


class translate : public Object
{
public:
	translate(Object* obj, const vec3& trans) :obj(obj), offset(trans) {}
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		ray trans_in = ray(in.start() - offset, in.direction());
		if (obj->hit(trans_in, min, max, param)) {
			param.pos = param.pos + offset;
			return true;
		}
		else return false;
	}
	vec3 center() const { return obj->center() + offset ; }
private:
	Object* obj;
	vec3 offset;
};


class rotate_y : public Object	// 绕中心点旋转
{
public:
	rotate_y(Object* obj, double angle) :obj(obj) {
		double radian = angle * M_PI / 180;
		sinTheta = sin(radian);
		cosTheta = cos(radian);
		obj_center = obj->center();
	}
	bool hit(const ray& in, double min, double max, hitParam& param) const {
		vec3 start = in.start() - obj_center;
		vec3 direction = in.direction();
		vec3 rotated_start = start, rotated_direction = direction;
		rotated_start[0] = cosTheta * start[0] - sinTheta * start[2];
		rotated_start[2] = sinTheta * start[0] + cosTheta * start[2];
		rotated_direction[0] = cosTheta * direction[0] - sinTheta * direction[2];
		rotated_direction[2] = sinTheta * direction[0] + cosTheta * direction[2];
		ray rotated_in = ray(rotated_start + obj_center, rotated_direction);
		if (obj->hit(rotated_in, min, max, param)) {
			vec3 pos = param.pos - obj_center;
			vec3 rotated_pos = pos;
			rotated_pos[0] = cosTheta * pos[0] + sinTheta * pos[2];
			rotated_pos[2] = -sinTheta * pos[0] + cosTheta * pos[2];
			param.pos = rotated_pos + obj_center;
			vec3 n = param.normal;
			param.normal[0] = cosTheta * n[0] + sinTheta * n[2];
			param.normal[2] = -sinTheta * n[0] + cosTheta * n[2];
			return true;
		}
		else return false;
	}
	vec3 center() const { return obj_center; }
private:
	Object* obj;
	vec3 obj_center;
	double sinTheta, cosTheta;
};

#endif