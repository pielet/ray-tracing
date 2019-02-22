#ifndef CAMERA
#define CAMERA

#include "vec3.h"
#include "ray.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

class Camera
{
public:
	Camera(vec3 pos, vec3 lookAt, vec3 up, double fov, double rate); // f = 1
	ray get_ray(double u, double v) {
		return ray(origin, left_bottom_point + horizontal * u + vertical * v - origin);
	}
private:
	vec3 origin;
	vec3 horizontal, vertical;
	vec3 left_bottom_point;
};

Camera::Camera(vec3 pos, vec3 lookAt, vec3 up, double fov, double rate) : origin(pos) {
	double theta = fov / 180 * M_PI;
	double half_height = tan(theta / 2);
	double half_width = half_height * rate;

	vec3 w = lookAt - pos;
	horizontal = normalize(w.cross(up)) * (2 * half_width);
	vertical = normalize(horizontal.cross(w)) * (2 * half_height);
	left_bottom_point = origin + normalize(w) - horizontal / 2 - vertical / 2;
}

#endif