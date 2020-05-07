#define _USE_MATH_DEFINES
#define GLM_ENABLE_EXPERIMENTAL
#define TASKING_TBB
//#define EMBREE_BACKFACE_CULLING 1

#include "ONB.h"

#include <cmath>

#include <embree3/rtcore.h>
#include <math.h>
#include <chrono>
#include <random>
#include <GL/glew.h>

#include <glm.hpp>
#include <gtx/normal.hpp>

#include <thread>
#include <mutex>
#include <iostream>
#include <map>

#include "ac3d.h"
#include "beckers.h"
#include "emissivity.h"

#include <windows.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <string>
#include <omp.h>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include <time.h>       /* time */
//MATLAB
#include <mat.h>
#include <matrix.h>
#include <mat.h>
#include <mex.h>   
#include <windows.h>
#include <process.h> 
#include "FreeImage.h"
#include <fstream>
#include <iostream>
typedef struct{
	float val;
	unsigned int gID;
	unsigned int pID;
} Result;

typedef struct{
	int r;
	int g;
	int b;
} colorInt;

enum REF { DIFFUSE, SPECULAR, GLOSSY };


#define STE_BOLZ 5.670367e-08
#define MAX_LEVELS 5
#define EPS 1e-4f
#define NUMTHREADS 8
#define SKY_SIZE 230 
#define NRAYS_GLOSSY 100
#define SPEC_N 500

int nRays = 100;
REF R = REF(2);

int hora;

RTCDevice device;
RTCScene eScene;
unsigned int triID;
unsigned int quadID;
int triOffset = 0;

std::vector<glm::vec3> sc_vertices;
std::vector<float> vertices_t;
std::vector<float> materials;

std::vector<std::vector<float>> tsky;


std::vector<int> sc_triangles;
std::vector<int> sc_quads;
//std::vector<std::tuple<GLuint, GLuint, GLfloat>> sparseFFs;
std::vector<glm::vec3> baris;
std::vector<glm::vec3> quadbaris;
std::vector<glm::vec3> tribaris;

std::vector<glm::vec3> normals;

std::mutex ffLock;

int N = 0;

int totalProcessed = 0;
int omp_get_thread_num();

/*
//Modest wood
float emisTable[91] = { 0, 0.06161, 0.12281, 0.18321, 0.24240, 0.30000, 0.35559, 0.40878, 0.45918, 0.50638,
0.55000, 0.58970, 0.62548, 0.65741, 0.68556, 0.71000, 0.73088, 0.74872, 0.76412, 0.77768,
0.79000, 0.80155, 0.81240, 0.82247, 0.83169, 0.84000, 0.84736, 0.85389, 0.85974, 0.86506,
0.87000, 0.87467, 0.87906, 0.88311, 0.88677, 0.89000, 0.89274, 0.89505, 0.89699, 0.89862,
0.90000, 0.90118, 0.90222, 0.90317, 0.90408, 0.90500, 0.90596, 0.90696, 0.90798, 0.90900,
0.91000, 0.91096, 0.91192, 0.91289, 0.91391, 0.91500, 0.91616, 0.91733, 0.91842, 0.91933,
0.92000, 0.92034, 0.92042, 0.92033, 0.92016, 0.92000, 0.91990, 0.91988, 0.91990, 0.91995,
0.92000, 0.92002, 0.92003, 0.92002, 0.92001, 0.92000, 0.91999, 0.91999, 0.91999, 0.91999,
0.92000, 0.92000, 0.92000, 0.92000, 0.92000, 0.92000, 0.91999, 0.91999, 0.91999, 0.91999, 0.92000 };*/
/*
//Elena concrete
float emisTable[91] = { 0, 0.118, 0.236, 0.354, 0.472, 0.590, 0.648, 0.706, 0.764, 0.822, 0.880, 0.888,
0.896, 0.904, 0.912, 0.920, 0.924, 0.928, 0.932, 0.936, 0.940, 0.944, 0.948,
0.952, 0.956, 0.960, 0.962, 0.964, 0.966, 0.968, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970 };
*/
//Experimental perspective bayona
/*float emisTable[91] = { 0, 0.40, 0.70, 0.85, 0.88, 0.90, 0.91, 0.915, 0.92, 0.93,
0.94, 0.942, 0.944, 0.946, 0.948, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95 };*/
float emisTable[91] = { 0, 0.400, 0.700, 0.820, 0.860, 0.880, 0.890, 0.905, 0.910, 0.915
, 0.920, 0.921, 0.922, 0.923, 0.924, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925
, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925, 0.925 };

float emisTableGlass[91] = { 0, 0.118, 0.236, 0.354, 0.472, 0.590, 0.648, 0.706, 0.764, 0.822, 0.880, 0.888,
0.896, 0.904, 0.912, 0.920, 0.924, 0.928, 0.932, 0.936, 0.940, 0.944, 0.948,
0.952, 0.956, 0.960, 0.962, 0.964, 0.966, 0.968, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970, 0.970,
0.970, 0.970 };

float emisTableWood[91] = { 0, 0.375, 0.675, 0.795, 0.835, 0.855, 0.865, 0.880,
0.885, 0.890, 0.895, 0.896, 0.897, 0.898, 0.899, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900 };


std::vector<glm::vec3> colormap;
//std::vector<float> tsky;


void intersectionFilter(const RTCFilterFunctionNArguments* args)
{
	/* avoid crashing when debug visualizations are used */
	if (args->context == nullptr) return;

	assert(args->N == 1);
	int* valid = args->valid;
	RTCRay* ray = (RTCRay*)args->ray;
	RTCHit* hit = (RTCHit*)args->hit;

	/* ignore inactive rays */
	if (valid[0] != -1) return;

	/* ignore hit if self hit */

	int globalID = hit->primID;
	if (hit->geomID == quadID)
		globalID += triOffset;

	int self = ray->flags;
	if (globalID == self)
		valid[0] = 0;
}



void loadSceneEmbree() {
	device = rtcNewDevice("threads=0");
	//	rtcSetDeviceProperty(device, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED, 1);
	//	ssize_t pr = rtcGetDeviceProperty(device, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED);
	//	std::cout << "BFC prop is " <<  pr<<  " \n";

	eScene = rtcNewScene(device);
	RTCBuffer vertices =
		rtcNewSharedBuffer(device, sc_vertices.data(),
		sizeof(glm::vec3) * sc_vertices.size());

	if (sc_triangles.size() > 0) {
		RTCGeometry triGeom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
		rtcSetGeometryBuffer(triGeom, RTC_BUFFER_TYPE_VERTEX, 0,
			RTC_FORMAT_FLOAT3, vertices, 0, sizeof(glm::vec3),
			sc_triangles.size());

		GLuint* indexT = (GLuint*)rtcSetNewGeometryBuffer(
			triGeom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			3 * sizeof(GLuint), sc_triangles.size() / 3);

		for (GLuint i = 0; i < sc_triangles.size(); i++) {
			indexT[i] = sc_triangles[i];
		}
		//filter to avoid self hit
		rtcSetGeometryIntersectFilterFunction(triGeom, intersectionFilter);

		rtcSetGeometryVertexAttributeCount(triGeom, 1);
		//		rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, hair_indices, 0, sizeof(unsigned int), 1);
		rtcSetSharedGeometryBuffer(triGeom, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT, vertices_t.data(), 0, sizeof(float), vertices_t.size());

		rtcCommitGeometry(triGeom);
		triID = rtcAttachGeometry(eScene, triGeom);
		triOffset = sc_triangles.size() / 3;
	}
	if (sc_quads.size() > 0) {
		RTCGeometry quadGeom =
			rtcNewGeometry(device, RTC_GEOMETRY_TYPE_QUAD);
		rtcSetGeometryBuffer(quadGeom, RTC_BUFFER_TYPE_VERTEX, 0,
			RTC_FORMAT_FLOAT3, vertices, 0, sizeof(glm::vec3),
			sc_quads.size());
		GLuint* indexQ = (GLuint*)rtcSetNewGeometryBuffer(
			quadGeom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4,
			4 * sizeof(GLuint), sc_quads.size() / 4);
		for (GLuint i = 0; i < sc_quads.size(); i++) {
			indexQ[i] = sc_quads[i];
		}
		//filter to avoid self hit
		rtcSetGeometryIntersectFilterFunction(quadGeom, intersectionFilter);

		rtcSetGeometryVertexAttributeCount(quadGeom, 1);
		//		rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, hair_indices, 0, sizeof(unsigned int), 1);
		rtcSetSharedGeometryBuffer(quadGeom, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT, vertices_t.data(), 0, sizeof(float), vertices_t.size());

		rtcCommitGeometry(quadGeom);
		quadID = rtcAttachGeometry(eScene, quadGeom);
	}

	rtcCommitScene(eScene);
}


void DeleteEverything() {
	sc_vertices.clear();
	sc_triangles.clear();
	sc_quads.clear();
	rtcReleaseScene(eScene);
	rtcReleaseDevice(device);
}

void load_scene(const char * file) {
	std::cout << "loading... ";
	if (ac3d::load_scene(file, sc_vertices, sc_triangles, sc_quads)<0) {
		std::cout << "error\n";
		return;
	}
	std::cout << "done\n";
	std::cout << "scene has " << sc_vertices.size() << " points\n";
	std::cout << "scene has " << sc_triangles.size() / 3 << " triangles\n";
	std::cout << "scene has " << sc_quads.size() / 4 << " quads\n";

	auto start = std::chrono::high_resolution_clock::now();
	loadSceneEmbree();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "time: " << elapsed.count() << "\n";
}

void saveImage(RGBQUAD *cs, int width, int height, std::string file) {
	FreeImage_Initialise();
	FIBITMAP* bitmap = FreeImage_Allocate(width, height, 24);
	if (bitmap == NULL)
		exit(1);
	for (int i = 0; i<width; i++)
		for (int j = 0; j<height; j++)
			//			FreeImage_SetPixelColor(bitmap, (width-1)-i, (height-1)-j, &cs[i*height + j]);
			FreeImage_SetPixelColor(bitmap, i, j, &cs[i*height + j]);

	if (!FreeImage_Save(FIF_PNG, bitmap, file.c_str(), 0)) // El ultimo parametro son flags. Dejar siempre en 0.
		exit(1);
	FreeImage_DeInitialise();
}

void loadMaterials(){
	std::cout << "loading materials...\n";
	std::ifstream ifs("../mat");
	float mat;
	materials.clear();
	while (ifs >> mat)
	{
		materials.push_back(mat);
	}
	std::cout << "total materials " << materials.size() << " ...\n";
}

void loadFieldByPoint(std::string file){
	std::cout << "loading field by points...\n";

	std::ifstream ifs(file);
	float t;

	vertices_t.clear();
	while (ifs >> t) {
		vertices_t.push_back(t);
	}
	std::cout << "listo\n";
}

void loadColormap(){
	std::ifstream ifs("../colormap");

	std::cout << "cargando colormap...\n";

	glm::vec3 c;

	colormap.clear();
	while (ifs >> c.x >> c.y >> c.z) {
		colormap.push_back(c);
	}
	std::cout << "listo\n";
}

/*void loadSkyTemps(){
std::stringstream ss;
ss << "../skyBeck/skyTemp" << (hora*6)+1 << ".csv";
std::string s = ss.str();
std::ifstream ifs(s);
std::cout << "cargando sky temps...\n";

float ts;

tsky.clear();
while (ifs >> ts) {
tsky.push_back(ts);
}
std::cout << "listo\n";

}*/
void loadSkyTemps(){
	std::stringstream ss;
	ss << "../tsky";
	std::string s = ss.str();
	std::ifstream ifs(s);
	std::cout << "loading sky temps...\n";

	tsky.resize(48);
	for (int i = 0; i < 48; i++)
		tsky[i].resize(10);

	for (int i = 0; i < 48; i++){
		for (int j = 0; j < 10; j++){
			ifs >> tsky[i][j];
		}
	}
	std::cout << "total sky temps " << tsky.size() << " ...\n";
}

int getColor(float t){
	float tt = (t - 273.15);

	if (tt < 10.0)
		tt = 10.0;
	else if (tt > 40.0)
		tt = 40.0;

	int ind = floor(((tt - 10) / 30)*(colormap.size() - 1));
	return ind;
}

int getColor2(float t){
	float tt = (t - 273.15);

	if (tt < -10.0)
		tt = -10.0;
	else if (tt > 30.0)
		tt = 30.0;

	int ind = floor(((tt + 10) / 40)*(colormap.size() - 1));
	return ind;
}


void saveData(float* data, int width, int height){
	std::ofstream ffs_file("../thermo/temps");
	for (int i = 0; i<width; ++i) {
		for (int j = 0; j<height; ++j)
			ffs_file << data[i*height + j] << ' ';
		ffs_file << '\n';
	}
	ffs_file.close();
}

Result getSpecularlyReflectedTemperature(glm::vec3 reflDir, glm::vec3 hitPoint, float globalID){
	Result r;

	RTCRay ray;
	ray.org_x = hitPoint.x; ray.org_y = hitPoint.y;	ray.org_z = hitPoint.z;
	ray.tfar = 1000.0f; ray.tnear = EPS; ray.time = 0; ray.id = 0;
	//avoid self hit
	ray.flags = globalID;
	ray.dir_x = reflDir.x; 	ray.dir_y = reflDir.y; ray.dir_z = reflDir.z;

	RTCRayHit query;
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	query.ray = ray;
	query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	query.hit.primID = RTC_INVALID_GEOMETRY_ID;
	rtcIntersect1(eScene, &context, &query);
	float t;
	if (query.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
		glm::vec3 color(0.0f, 0.0f, 0.0f);

		r.gID = query.hit.geomID;
		r.pID = query.hit.primID;
		float u = query.hit.u;
		float v = query.hit.v;

		{
			rtcInterpolate0(rtcGetGeometry(eScene, r.gID), r.pID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &t, 1);
		}
	}
	else {
		float angle = asin(reflDir.z) * 180 / M_PI;
		if (angle > 0 && angle < 90){
			int angleL = floor(angle / 10);
			int angleU = ceil(angle / 10);
			float rest = angle / 10 - angleL;
			t = (1 - rest)*tsky[hora * 2][9 - angleL] + rest*tsky[hora * 2][9 - angleU];//- 5.
		}
		else
			t = 275.0;

		/*
		int skyTile = beckers(reflDir);
		if (skyTile != 0)
		t = tsky[skyTile - 1];
		else{
		t = 283.0;
		printf("ERROR \n");
		}*/
	}
	r.val = t;
	return r;
}

float getDiffuselyReflectedTemperature(glm::vec3 orig, ONB base, float rings, float	Drho){
	float tilespering = 1;
	float reflT = 0;
	float t2;
	for (float ring = 1; ring <= 2 * rings + 1; ring += 2) {
		float Dtheta = 2 * M_PI / tilespering;
		float rho = (ring - 1.0)*Drho;
		for (float k = 1; k <= tilespering; k++) {
			float theta = (k - 0.5)*Dtheta;
			float x = rho* cos(theta);
			float y = rho* sin(theta);
			float z = sqrt(1.0 - x*x - y*y);

			RTCRay ray;
			ray.org_x = orig.x; ray.org_y = orig.y; ray.org_z = orig.z;
			ray.tfar = 1000.0f; ray.tnear = 0; ray.time = 0; ray.flags = 0;
			ray.id = 0;
			ray.flags = -1; //pass id of origin to filter
			RTCRayHit query;
			RTCIntersectContext context;

			glm::vec3 dir = base.LocalToWorld(glm::vec3(x, y, z));
			dir = glm::normalize(dir);
			ray.dir_x = dir.x; ray.dir_y = dir.y; ray.dir_z = dir.z;
			rtcInitIntersectContext(&context);

			query.ray = ray;
			query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
			query.hit.primID = RTC_INVALID_GEOMETRY_ID;

			rtcIntersect1(eScene, &context, &query);
			if (query.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
				glm::vec3 color(0.0f, 0.0f, 0.0f);

				unsigned int gID = query.hit.geomID;
				unsigned int pID = query.hit.primID;
				float u = query.hit.u;
				float v = query.hit.v;
				{
					rtcInterpolate0(rtcGetGeometry(eScene, gID), pID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &t2, 1);
				}
				reflT += t2;
			}
			else{
				float t;
				float angle = asin(dir.z) * 180 / M_PI;
				if (angle > 0 && angle < 90){
					int angleL = floor(angle / 10);
					int angleU = ceil(angle / 10);
					float rest = angle / 10 - angleL;
					t = (1 - rest)*tsky[hora * 2][9 - angleL] + rest*tsky[hora * 2][9 - angleU];
				}
				else
					t = 283.0;
				reflT += t;
				/*int skyTile = beckers(dir);
				if (skyTile != 0)
				reflT += tsky[skyTile - 1];
				else
				reflT += 283.0;*/
			}
		}
		tilespering = (ring + 1) * 4;
	}
	return reflT / nRays;
}

float cosAngleBetween(glm::vec3 v1, glm::vec3 v2){
	float dot = glm::dot(v1, v2);
	float lenSq1 = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
	float lenSq2 = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
	return (dot / sqrt(lenSq1 * lenSq2));
	//	return acos(dot / sqrt(lenSq1 * lenSq2));
}

//http://www.cim.mcgill.ca/~derek/ecse689_a3.html
vec3 specLobe(ONB base, float e1, float e2, float N){
	float pp = pow(e1, (2 / (N + 1)));
	float x = sqrt(1 - pp)*cos(2 * M_PI*e2);
	float y = sqrt(1 - pp)*sin(2 * M_PI*e2);
	float z = pow(e1, (1 / (N + 1)));

	vec3 vec(x, y, z);

	vec = base.LocalToWorld(normalize(vec));
	return normalize(vec);;
}

float RadicalInverse_VdC(uint bits)
{
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

vec2 Hammersley(uint i, uint N)
{
	return vec2(float(i) / float(N), RadicalInverse_VdC(i));
}

#include <chrono>  // for high_resolution_clock

void generateThermography(){
	std::stringstream ss;
	ss << "../chp_T/T" << hora;
	std::string s = ss.str();
	loadFieldByPoint(s);
	load_scene("C:\\Users\\Jose\\Desktop\\Setiembre19\\Embree\\embree-rendering\\bayonneBox.ac");
	
	loadSkyTemps();

	//CAMERA DEFINTION
	//	glm::vec3 camOrig = glm::vec3(9.04,86,0.7);
	//	glm::vec3 camDir = glm::vec3(-tan(11 * M_PI / 180), -1, tan(11 * M_PI / 180));
	//glm::vec3 camOrig = glm::vec3(7.87, 86, 1.5);
	//glm::vec3 camDir = glm::vec3(-tan(8.3 * M_PI / 180), -1, tan(9.8 * M_PI / 180));
	//camOrig = camOrig - camDir*1.2f;

	glm::vec3 recta = glm::vec3(-tan(11 * M_PI / 180), -1, 0);


	//3 , 10
	glm::vec3 camOrig = glm::vec3(6.92,77.3,1);
	glm::vec3 camDir = glm::vec3(-tan(11 * M_PI / 180), -1, tan(11.3 * M_PI / 180));
	camDir = glm::normalize(camDir);

	camOrig = camOrig - recta*8.0f;

	//	glm::vec3 camOrig = glm::vec3(-16.2,34.4,2.52);
	//	glm::vec3 camDir = glm::vec3(1,0,0);

	glm::vec3 up = glm::vec3(0, 0, 1);
	glm::vec3 right = glm::cross(camDir, up);
	right = glm::normalize(right);
	float fovU = (25)*M_PI / 180;
	float fovR = (18.7)*M_PI / 180;
	float aspectRatio = fovU / fovR;
	std::cout << aspectRatio << "\n";

	//IMAGE DEFINITION
	int width = 180;
	const int AA = 4; //4x4 window anti-aliasing
	//int width = 500;
	//int width = 180*4;
	//int width = 150;
	int height = width*aspectRatio;
	RGBQUAD *apColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *reColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *diffColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *emisColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	float *tempData = (float*)malloc(width*height*sizeof(float));


	int countMAL = 0;
	int countBIEN = 0;

	//Beckers subdivision
	float rings = ceil((-1 + sqrt(nRays)) / 2);
	nRays = 1 + 4 * rings*(rings + 1); // n >= al n original
	float	Drho = 1 / (2 * rings + 1); //Delta rho
	//
	totalProcessed = 0;
	//COMPUTE
	omp_set_num_threads(NUMTHREADS);
	float* reflTAnt = (float*)malloc(NUMTHREADS*sizeof(float));

	colorInt   apAA[NUMTHREADS][AA*AA];
	colorInt   reAA[NUMTHREADS][AA*AA];
	colorInt diffAA[NUMTHREADS][AA*AA];
	colorInt emisAA[NUMTHREADS][AA*AA];

	float dfovR = fovR/width;
	float dfovU = fovR/height;

#pragma omp parallel for 
	for (int i = 0; i < width; i++){
		ffLock.lock();
		if (totalProcessed % (20) == 0)
			std::cout << 100.0f*totalProcessed / width << "%\n";
		totalProcessed++;
		ffLock.unlock();

		for (int j = 0; j < height; j++){
			float normalized_i = ((float)i / width) - 0.5;
			float normalized_j = ((float)j / height) - 0.5;

			for (int iAA = 0; iAA < AA; iAA++){
				for (int jAA = 0; jAA < AA; jAA++){
					float normalized_iAA = ((float)iAA / AA) - 0.5;
					float normalized_jAA = ((float)jAA / AA) - 0.5;

					//first, move to pixel center in camera space
					glm::vec3 pixelLocation = normalized_i * 2 * (right*tan(fovR / 2)) + normalized_j * 2 * (up*tan(fovU / 2));
					//second, move to cuadrant inside pixel due to AA
					pixelLocation += normalized_iAA * 2 * (right*tan(dfovR / 2)) + normalized_jAA * 2 * (up*tan(dfovU / 2));
					//translate to cam pos and dir
					pixelLocation += camOrig + camDir;

					glm::vec3 dir = pixelLocation - camOrig;
					dir = glm::normalize(dir);

					RTCRay ray;
					ray.org_x = camOrig.x; ray.org_y = camOrig.y;	ray.org_z = camOrig.z;
					ray.tfar = 1000.0f; ray.tnear = EPS; ray.time = 0; ray.id = 0; ray.flags = -1;
					ray.dir_x = dir.x; 	ray.dir_y = dir.y; ray.dir_z = dir.z;

					RTCRayHit query;
					RTCIntersectContext context;
					rtcInitIntersectContext(&context);
					query.ray = ray;
					query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
					query.hit.primID = RTC_INVALID_GEOMETRY_ID;
					rtcIntersect1(eScene, &context, &query);

					float t = 0;
					float aparentT = 0;
					if (query.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
						int globalID = query.hit.primID;
						if (query.hit.geomID == quadID)
							globalID += triOffset;

						glm::vec3 apColor(0.0f, 0.0f, 0.0f);
						glm::vec3 reColor(0.0f, 0.0f, 0.0f);
						glm::vec3 refColor(0.0f, 0.0f, 0.0f);

						unsigned int gID = query.hit.geomID;
						unsigned int pID = query.hit.primID;
						float u = query.hit.u;
						float v = query.hit.v;

						{
							rtcInterpolate0(rtcGetGeometry(eScene, gID), pID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &t, 1);
						}

						glm::vec3 originalDir = glm::vec3(query.ray.dir_x, query.ray.dir_y, query.ray.dir_z);
						glm::vec3 hitNormal = glm::normalize(glm::vec3(query.hit.Ng_x, query.hit.Ng_y, query.hit.Ng_z));
						dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
						glm::vec3 orig = glm::vec3(query.ray.org_x, query.ray.org_y, query.ray.org_z) + query.ray.tfar * originalDir + hitNormal*EPS;
						float specN = SPEC_N;

						//float * emiT = emisTable;
						//if (materials[globalID] != 3.0)
						//	emiT = emisTableWood;

						float * emiT = emisTable;
						if (materials[globalID] == 2.0 || materials[globalID] == 8.0) // wood and window
						{
							emiT = emis10;
							specN = 200;
						}
						else if (materials[globalID] == 3.0) //wall
						{
							emiT = emis15;
							specN = 300;
						}
						else if (materials[globalID] == 5.0) //roof
						{
							emiT = emis50;
							specN = 1000;
						}
						else if (materials[globalID] == 7.0) //street
						{
							emiT = emis30;
							specN = 400;
						}
						else if (materials[globalID] == 4.0) //stone
						{
							emiT = emis10;
							specN = 500;
						}

						if (orig.y > 45. && orig.x<0.173 && orig.x>-0.951 && (materials[globalID] == 8.0)){
							emiT = emisTableGlass;
							specN = 1e5;
						}


						ONB base(hitNormal);
						glm::vec3 dirLocal = base.WorldToLocal(-originalDir);

						dirLocal = glm::normalize(dirLocal);
						float angulo = (asin(dirLocal.z) * 180 / M_PI);

						if (angulo < 0){
							apColor = colormap[getColor(t)];
							countMAL++;
						}
						else{
							float reflT = 0;
							if (R == DIFFUSE)
								reflT = getDiffuselyReflectedTemperature(orig, base, rings, Drho);
							else if (R == SPECULAR){
								dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
								unsigned int gID_hitted;
								unsigned int pID_hitted;
								Result res;
								res = getSpecularlyReflectedTemperature(dir, orig, globalID);
								reflT = res.val;
							}
							else if (R == GLOSSY){
								dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
								float theta = (angulo < 5 ? angulo : 5) * M_PI / 180.0;
								float sinAngleMAX = sin(theta);
								int nRaysX = round(sqrt(NRAYS_GLOSSY));
								int count = 0;
								ONB base(normalize(dir));
								for (int rx = 0; rx < nRaysX; rx++)
									for (int ry = 0; ry < nRaysX; ry++){
										//float x = (float)(rx+1) / nRaysX;
										//float y = (float)(ry+1) / nRaysX;
										//float x = (float)rand() / RAND_MAX;
										//float y = (float)rand() / RAND_MAX;
										vec2 Xi = Hammersley(rx*nRaysX + ry, nRaysX*nRaysX);
										float x = Xi.x;
										float y = Xi.y;

										vec3 perturbed = specLobe(base, x, y, specN);
										Result res = getSpecularlyReflectedTemperature(perturbed, orig, -1);
										if ((res.gID == gID) && (res.pID == pID)){
											reflT += res.val;
											count++;
										}
										else{
											reflT += res.val;
											count++;
										}

									}
								if (count > 0){
									reflT = reflT / (count);
									reflTAnt[omp_get_num_threads()] = reflT;
								}
								else{
									reflT = reflTAnt[omp_get_num_threads()];
								}
							}
							int ang1 = floor(angulo);
							int ang2 = ceil(angulo);
							float coef = angulo - ang1;

							float emis = (1 - coef)*emiT[ang1] + coef*emiT[ang2];
							//emis = 0.95;
							float refl = 1 - emis;
							float directFlux = emis*pow(t, 4);
							float reflectedFlux = refl*pow(reflT, 4.0);
							aparentT = pow(directFlux + reflectedFlux, 1.0 / 4.0);

							apColor = colormap[getColor(aparentT)];
							refColor = colormap[getColor2(reflT)];

							countBIEN++;

							tempData[i*height + j] = aparentT;

							emisAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(emis * 255);
							emisAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(emis * 255);
							emisAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(emis * 255);
						}
						apAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(apColor.r * 255);
						apAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(apColor.g * 255);
						apAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(apColor.b * 255);

						reColor = colormap[getColor(t)];
						reAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(reColor.r * 255);
						reAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(reColor.g * 255);
						reAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(reColor.b * 255);

						diffAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(refColor.r * 255);
						diffAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(refColor.g * 255);
						diffAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(refColor.b * 255);

					}
					else{
						apAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(255);
						apAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(255);
						apAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(255);

						reAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(255);
						reAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(255);
						reAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(255);

						diffAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(255);
						diffAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(255);
						diffAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(255);

						emisAA[omp_get_thread_num()][iAA*AA + jAA].r = (int)(255);
						emisAA[omp_get_thread_num()][iAA*AA + jAA].g = (int)(255);
						emisAA[omp_get_thread_num()][iAA*AA + jAA].b = (int)(255);
					}


				}
			}
			//gather Anti Aliasing
			colorInt ap = { 0, 0, 0 }, re = { 0, 0, 0 }, diff = { 0, 0, 0 }, emis = { 0, 0, 0 };
			for (int iAA = 0; iAA < AA; iAA++){
				for (int jAA = 0; jAA < AA; jAA++){
					ap.r += apAA[omp_get_thread_num()][iAA*AA + jAA].r;
					ap.g += apAA[omp_get_thread_num()][iAA*AA + jAA].g;
					ap.b += apAA[omp_get_thread_num()][iAA*AA + jAA].b;

					re.r += reAA[omp_get_thread_num()][iAA*AA + jAA].r;
					re.g += reAA[omp_get_thread_num()][iAA*AA + jAA].g;
					re.b += reAA[omp_get_thread_num()][iAA*AA + jAA].b;

					diff.r += diffAA[omp_get_thread_num()][iAA*AA + jAA].r;
					diff.g += diffAA[omp_get_thread_num()][iAA*AA + jAA].g;
					diff.b += diffAA[omp_get_thread_num()][iAA*AA + jAA].b;

					emis.r += emisAA[omp_get_thread_num()][iAA*AA + jAA].r;
					emis.g += emisAA[omp_get_thread_num()][iAA*AA + jAA].g;
					emis.b += emisAA[omp_get_thread_num()][iAA*AA + jAA].b;
				}
			}

			apColors[i*height + j].rgbRed    =  (BYTE) (ap.r / (AA*AA));
			apColors[i*height + j].rgbGreen  =  (BYTE) (ap.g / (AA*AA));
			apColors[i*height + j].rgbBlue   =  (BYTE) (ap.b / (AA*AA));
											   
			reColors[i*height + j].rgbRed    =  (BYTE) (re.r / (AA*AA));
			reColors[i*height + j].rgbGreen  =  (BYTE) (re.g / (AA*AA));
			reColors[i*height + j].rgbBlue   =  (BYTE) (re.b / (AA*AA));
											   
			diffColors[i*height + j].rgbRed  =  (BYTE) (diff.r / (AA*AA));
			diffColors[i*height + j].rgbGreen=  (BYTE) (diff.g / (AA*AA));
			diffColors[i*height + j].rgbBlue =  (BYTE) (diff.b / (AA*AA));
											   
			emisColors[i*height + j].rgbRed  =  (BYTE) (emis.r / (AA*AA));
			emisColors[i*height + j].rgbGreen=  (BYTE) (emis.g / (AA*AA));
			emisColors[i*height + j].rgbBlue =  (BYTE) (emis.b / (AA*AA));
			
		}
	}

	if (hora == 14 || hora == 23 )
		saveData(tempData, width, height);

	std::cout << "countBIEN, countMAL, total " << countBIEN << " " << countMAL << " " << width*height << "\n";
	//SAVE IMAGE
	std::stringstream ss1;
	ss1 << "../thermo/real" << hora << ".png";

	std::stringstream ss2;
	ss2 << "../thermo/apparent" << hora << ".png";

	std::stringstream ss3;
	ss3 << "../thermo/refl" << hora << ".png";

	std::stringstream ss4;
	ss4 << "../thermo/emis.png";

	saveImage(reColors, width, height, ss1.str());
	saveImage(apColors, width, height, ss2.str());
	saveImage(diffColors, width, height, ss3.str());

	if (hora == 0)
		saveImage(emisColors, width, height, ss4.str());

	free(reColors);
	free(apColors);
	free(diffColors);
};

int main() {
	srand(static_cast <unsigned> (time(0)));

	loadColormap();
	loadMaterials();
	std::cout << "size colormap " << colormap.size() << "\n";

	for (int i = 0; i < 24; i++){
		hora = i;
		std::cout << "Generating hour " << hora << " ...\n";
		generateThermography();
		DeleteEverything();
	}


	//std::cin >> nRays;
	return 0;
}