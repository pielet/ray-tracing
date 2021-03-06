#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "vec3.h"
#include "ray.h"
#include "camera.h"
#include "object.h"
#include "material.h"

using namespace std;

vec3 color(const ray& in, Object_list object, vec3 ambient, int depth);
vec3 Phong_shading(const ray& in, Object_list objects, double ambient_rate, vec3 light_pos);

int main()
{
	int nx = 1200;
	int ny = 800;

	Object_list obj_list;
	
	// scene1
	Camera cam(vec3(-0.3, 0.6, 0.5), vec3(0, 0.4, 0), vec3(0, 1, 0), 80, double(nx)/ ny);
	obj_list.add(new Sphere(vec3(0.3, -0.3, -0.3), 0.2, new transparent(1.5), new constant(vec3(1, 1, 1))));
	obj_list.add(new Sphere(vec3(0, 0, -1), 0.5,        new fuzzy(),          new image("texture\\sphere.jpg")));
	// obj_list.add(new Sphere(vec3(0, -100.5, -1), 100, new fuzzy(),         new constant(vec3(0.5, 0.2, 0.8))));
	obj_list.add(new Sphere(vec3(0, -100.5, -1), 100,   new fuzzy(),          new constant(vec3(0, 0.75, 1))));
	obj_list.add(new Sphere(vec3(1, 0.6, -0.6), 0.3,    new specular(0),      new constant(vec3(0.8, 0.8, 0.8))));
	obj_list.add(new cube(vec3(0.8, -0.5, -1), 0.8,     new fuzzy(),          new image("texture\\rect.jpg")));
	
	obj_list.add(new rect(1, 3, -2, 2, -2, 0, new light(), new constant(vec3(1, 1, 1))));
	
	/*
	// scene2
	Camera cam(vec3(-0.1, 0.7, 1.2), vec3(0, 0.3, 0), vec3(0, 1, 0), 60, double(nx)/ ny);
	cout << "start loading obj file ...\n";
	obj_list.add(new model("model\\beautiful bunny\\bunny.obj", 4, new fuzzy(), new image("model\\beautiful bunny\\white.jpg")));
	obj_list.add(new rect(1, 0, -5, 5, -3, 5, new fuzzy(), new chessboard(vec3(0.5, 0.8, 0.1), 30)));
	cout << "finish\n";
	*/
	/*
	// scene 3
	Camera cam(vec3(0, 0, 1), vec3(0, 0, 0), vec3(0, 1, 0), 60, double(nx) / ny);

	obj_list.add(new Sphere(vec3(0, 0, -1), 0.5, new fuzzy(), new constant(vec3(0, 0.75, 1))));
	obj_list.add(new Sphere(vec3(0, -100.5, -1), 100, new fuzzy(), new constant(vec3(0.75, 0.75, 0))));
	*/

	vec3 ambient(0, 0, 0);
	int sample_n = 8;

	srand(unsigned(time(nullptr)));
	ofstream f("image\\test.ppm");
	if (f.is_open()) { 
		f << "P3\n" << nx << " " << ny << " 255\n";
		double u, v;
		for (int i = ny - 1; i >= 0; --i) {
			cout << i << endl;
			for (int j = 0; j < nx; ++j) {
				vec3 pix_color(0, 0, 0);
				for (int n = 0; n < sample_n; ++n) {
					u = (j + double(rand()) / RAND_MAX) / nx;
					v = (i + double(rand()) / RAND_MAX) / ny;
					ray in = cam.get_ray(u, v);
					pix_color = pix_color + color(in, obj_list, ambient, 50);
					// pix_color = pix_color + Phong_shading(in, obj_list, 0.05, vec3(-5, 5, 0));
				}
				pix_color = pix_color / sample_n;
				pix_color = vec3(sqrt(pix_color[0]), sqrt(pix_color[1]), sqrt(pix_color[2]));
				f << pix_color << endl;
			}
			f << endl;
		}
		f.close();
	}

	return 0;
}


vec3 color(const ray& in, Object_list objects, vec3 ambient, int depth) {
	if (depth == 0) return vec3(0, 0, 0);

	hitParam param;
	int _;
	if (objects.hit(in, 0.001, DBL_MAX, param, _)) {
		vec3 obj_color = param.tp->value(param.u, param.v);
		if (param.mp->is_light())
			return obj_color;
		else {
			ray out;
			if (param.mp->scatter(in, param, out))
				//return color(out, objects, ambient, depth - 1) * 0.5;
				return obj_color * color(out, objects, ambient, depth - 1);
			else return vec3(0, 0, 0);
		}
	}
	return ambient;
	
	// vec3 unit_direction = normalize(in.direction());
	// float t = 0.5*(unit_direction[1] + 1.0);
	// return vec3(1.0, 1.0, 1.0)*(1-t) + vec3(0.5, 0.7, 1.0)*t;
	
}


vec3 Phong_shading(const ray& in, Object_list objects, double ambient_rate, vec3 light_pos) {
	hitParam param, tmp_p;
	int _;
	if (objects.hit(in, 0.001, DBL_MAX, param, _)) {
		vec3 obj_color = param.tp->value(param.u, param.v);
		vec3 hit_point_2_light = light_pos - param.pos;

		// only ambient
		if (objects.hit(ray(param.pos, hit_point_2_light), 0.001, DBL_MAX, tmp_p, _))
			return obj_color * ambient_rate;

		// diffuse
		double cos = param.normal.dot(hit_point_2_light) / param.pos.length() / hit_point_2_light.length();
		if (cos < 0.01) cos = 0;
		vec3 diffuse_color = obj_color * cos;
		// specular
		ray reflect;
		vec3 specular_color;
		if (specular(0).scatter(in, param, reflect)) {
			cos = reflect.direction().dot(hit_point_2_light) / reflect.direction().length() / hit_point_2_light.length();
			cos = pow(cos, 5);
			specular_color = vec3(cos, cos, cos);
		}

		vec3 total_color = obj_color * ambient_rate + diffuse_color + specular_color;

		for (int i = 0; i < 3; ++i) {
			if (total_color[i] > 1) total_color[i] = 1;
		}

		return total_color;
	}
	else return vec3(0, 0, 0);
}


/*
vec3 color(const ray& in, Object_list objects, vec3 ambient, int depth, vec3 light_pos) {
if (depth == 0) return vec3(0, 0, 0);

hitParam param, tmp_p;
int _;
if (objects.hit(in, 0.001, DBL_MAX, param, _)) {
vec3 obj_color = param.tp->value(param.u, param.v);

// Phong
vec3 hit_point_2_light = light_pos - param.pos;
vec3 ambient_color = ambient_color = obj_color * 0.05;
vec3 diffuse_color, specular_color;

if (!objects.hit(ray(param.pos, hit_point_2_light), 0.001, DBL_MAX, tmp_p, _)) {
// diffuse
double cos = param.normal.dot(hit_point_2_light) / param.pos.length() / hit_point_2_light.length();
if (cos < 0.01) cos = 0;
diffuse_color = obj_color * cos;
// specular
ray reflect;
if (specular(0).scatter(in, param, reflect)) {
cos = reflect.direction().dot(hit_point_2_light) / reflect.direction().length() / hit_point_2_light.length();
cos = pow(cos, 5);
if (cos < 0.01) cos = 0;
specular_color = vec3(cos, cos, cos);
}
}
vec3 total_color = ambient_color + diffuse_color + specular_color;

ray out;
if (param.mp->scatter(in, param, out)) {
vec3 return_color = color(out, objects, ambient, depth - 1, light_pos)*0.5;
if(return_color[0] > 0 || return_color[1] > 0 || return_color[2] > 0)
total_color = total_color * color(out, objects, ambient, depth - 1, light_pos);
}

for (int i = 0; i < 3; ++i) {
if (total_color[i] > 1) total_color[i] = 1;
}
return total_color;
}
return ambient;
}
*/