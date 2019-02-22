// 基础类，定义了三维向量
#ifndef VEC3
#define VEC3
#include <iostream>
#include <cmath>
using namespace std;

class vec3
{
public:
	vec3(double x, double y, double z) { pos[0] = x; pos[1] = y; pos[2] = z; }
	vec3(const vec3& v) { pos[0] = v[0]; pos[1] = v[1]; pos[2] = v[2]; }
	vec3() { pos[0] = pos[1] = pos[2] = 0; }
	// 输入输出
	friend istream& operator>>(istream& in, vec3& t);
	friend ostream& operator<<(ostream& out, const vec3& t);
	// 下标和向量间运算
	double& operator[](unsigned i) { return pos[i]; }
	double operator[](unsigned i) const { return pos[i]; }
	vec3 operator- () const { return vec3(-pos[0], -pos[1], -pos[2]); }
	vec3 operator+ (const vec3 &v) const { return vec3(v.pos[0] + pos[0], v.pos[1] + pos[1], v.pos[2] + pos[2]); }
	vec3 operator- (const vec3& v) const { return vec3(pos[0]- v.pos[0], pos[1] - v.pos[1], pos[2] - v.pos[2]); }
	vec3 operator* (const vec3& v) const { return vec3(v.pos[0] * pos[0], v.pos[1] * pos[1], v.pos[2] * pos[2]); }
	vec3& operator= (const vec3& v) {
		for (int i = 0; i < 3; ++i)
			pos[i] = v.pos[i];
		return *this;
	}
	// 和标量运算
	vec3 operator* (double k) const { return vec3(pos[0] * k, pos[1] * k, pos[2] * k); }
	vec3 operator/ (double k) const { return vec3(pos[0] / k, pos[1] / k, pos[2] / k); }
	// 点乘/叉乘
	double dot(const vec3& v) const { return v.pos[0] * pos[0] + v.pos[1] * pos[1] + v.pos[2] * pos[2]; }
	vec3 cross(const vec3& v) const {
		return vec3(pos[1] * v.pos[2] - pos[2] * v.pos[1],
			pos[2] * v.pos[0] - pos[0] * v.pos[2],
			pos[0] * v.pos[1] - pos[1] * v.pos[0]);
	}
	// 长度
	double squared_length() const { return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]; }
	double length() const { return sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]); }
	void normalize() {
		double k = 1 / this->length();
		for (int i = 0; i < 3; ++i)
			pos[i] *= k;
	}
	// 最值
	double max() const {
		double a = ((pos[0] > pos[1]) ? pos[0] : pos[1]);
		return ((a > pos[2]) ? a : pos[2]);
	}
private:
	double pos[3];
};

inline istream& operator>>(istream& in, vec3& t) {
	in >> t.pos[0] >> t.pos[1] >> t.pos[2];
	return in;
}

inline ostream& operator<<(ostream& out, const vec3& t) {
	for (int i = 0; i < 3; ++i)
		out << int(255.99*t[i]) << " ";
	return out;
}

inline vec3 normalize(const vec3& t) {
	return t / t.length();
}

inline double distance(const vec3& v1, const vec3& v2) {
	return sqrt((v1[0] - v2[0])*(v1[0] - v2[0]) +
		(v1[1] - v2[1])*(v1[1] - v2[1]) +
		(v1[2] - v2[2])*(v1[2] - v2[2]));
}
 
#endif // !VEC3
