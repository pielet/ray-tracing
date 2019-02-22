#ifndef RAY
#define RAY

#include "vec3.h"

class ray
{
public:
	ray(const vec3& start, const vec3& direction) : A(start), B(direction) {};
	ray() {}
	vec3 position(double t) const { return A + B * t; }
	vec3 start() const { return A; }
	vec3 direction() const { return B; }
private:
	vec3 A, B;
};

#endif