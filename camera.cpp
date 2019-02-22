#include "camera.h"

Camera::Camera(vec3 pos, vec3 lookAt, vec3 up, double fov, double rate) : origin(pos) {
	double theta = fov / 180 * PI;
	double half_height = tan(theta / 2);
	double half_width = half_height * rate;

	vec3 w = lookAt - pos;
	horizontal = normalize(w.cross(up)) * (2 * half_width);
	vertical = normalize(horizontal.cross(w)) * (2 * half_height);
	left_bottom_point = origin + normalize(w) - horizontal / 2 - vertical / 2;
}

ray Camera::get_ray(double u, double v) {
	return ray(origin, left_bottom_point + horizontal * u + vertical * v - origin);
}