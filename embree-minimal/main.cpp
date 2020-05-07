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
#include "UCDimporter.h"
#include "colormap.h"


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
#include "EmbreeUtils.h"
#include "settingsLoader.h"
#include "materialLoader.h"

#include <fstream>
#include <iostream>

const float STE_BOLZ = 5.670367e-08;
const float EPS = 1e-4f;
const int NUMTHREADS = 8;

int NRAYS_GLOSSY;

std::mutex ffLock;
int totalProcessed = 0;
int omp_get_thread_num();

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

void saveData(float* data, int width, int height){
	std::ofstream ffs_file("../results/temps");
	for (int i = 0; i<width; ++i) {
		for (int j = 0; j<height; ++j)
			ffs_file << data[i*height + j] << ' ';
		ffs_file << '\n';
	}
	ffs_file.close();
}

Result getSpecularlyReflectedTemperature(float* tsky, glm::vec3 reflDir, glm::vec3 hitPoint, float globalID){
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
			t = (1 - rest)*tsky[9 - angleL] + rest*tsky[9 - angleU];//- 5.
		}
		else
			t = 275.0;
	}
	r.val = t;
	return r;
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


void printProgress(double percentage) {
	if (percentage > 0.95) percentage=1.0;
	int barLength = 60;	float pos = percentage * barLength;
	std::cout << "[";
	for (int i = 0; i != barLength; ++i){
		if (i < pos) std::cout << "|";
		else std::cout << " "; 	}
	std::cout << "] " << round(percentage*100) <<"% \r";
	if (percentage > 0.95) std::cout << "\n";
}

void generateThermography(float*tsky, std::vector<int> &matIDs, settings &s, material *matProps){
	std::cout << "Rendering thermography... \n";
	glm::vec3 camOrig = s.cameraCenter;
	glm::vec3 camDir = s.cameraDirection;
	camDir = glm::normalize(camDir);
	glm::vec3 up = s.cameraUp;
	glm::vec3 right = glm::cross(camDir, up);
	right = glm::normalize(right);
	float fovU = s.fovVertical*(M_PI / 180.0f);
	int width = s.imageWidth;
	int height = s.imageHeight;
	float aspectRatio = (float)s.imageHeight / (float)s.imageWidth;

	float fovR = fovU / aspectRatio;
	const int AA = round(sqrt(s.aa));

	RGBQUAD *apColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *reColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *diffColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *emisColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	float *tempData = (float*)malloc(width*height*sizeof(float));

	int countMAL = 0;
	int countBIEN = 0;

	totalProcessed = 0;
	//COMPUTE
	omp_set_num_threads(NUMTHREADS);
	float* reflTAnt = (float*)malloc(NUMTHREADS*sizeof(float));

	colorInt* apAA = (colorInt*)malloc(sizeof(colorInt)*NUMTHREADS*AA*AA);
	colorInt* reAA = (colorInt*)malloc(sizeof(colorInt)*NUMTHREADS*AA*AA);
	colorInt* diffAA = (colorInt*)malloc(sizeof(colorInt)*NUMTHREADS*AA*AA);
	colorInt* emisAA = (colorInt*)malloc(sizeof(colorInt)*NUMTHREADS*AA*AA);

	float dfovR = fovR/width;
	float dfovU = fovR/height;
#pragma omp parallel for 
	for (int i = 0; i < width; i++){
		ffLock.lock();
		if (totalProcessed % (width/30) == 0)
			printProgress(totalProcessed / (float)width);
			//std::cout << 100.0f*totalProcessed / width << "%\n";
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

						float * emiT = matProps[matIDs[globalID]].emisTable;
						float specN = matProps[matIDs[globalID]].specular_lobe_size;

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
					
							dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
							float theta = (angulo < 5 ? angulo : 5) * M_PI / 180.0;
							float sinAngleMAX = sin(theta);
							int nRaysX = round(sqrt(NRAYS_GLOSSY));
							int count = 0;
							ONB base(normalize(dir));
							for (int rx = 0; rx < nRaysX; rx++)
								for (int ry = 0; ry < nRaysX; ry++){
									vec2 Xi = Hammersley(rx*nRaysX + ry, nRaysX*nRaysX);
									float x = Xi.x;
									float y = Xi.y;

									vec3 perturbed = specLobe(base, x, y, specN);
									Result res = getSpecularlyReflectedTemperature(tsky,perturbed, orig, -1);
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

							emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r = (int)(emis * 255);
							emisAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].g = (int)(emis * 255);
							emisAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].b = (int)(emis * 255);
						}
						apAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].r = (int)(apColor.r * 255);
						apAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].g = (int)(apColor.g * 255);
						apAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].b = (int)(apColor.b * 255);

						reColor = colormap[getColor(t)];
						reAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].r = (int)(reColor.r * 255);
						reAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].g = (int)(reColor.g * 255);
						reAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].b = (int)(reColor.b * 255);

						diffAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].r = (int)(refColor.r * 255);
						diffAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].g = (int)(refColor.g * 255);
						diffAA[omp_get_thread_num()*AA*AA + iAA*AA + jAA].b = (int)(refColor.b * 255);

					}
					else{
						apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r = (int)(255);
						apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g = (int)(255);
						apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b = (int)(255);

						reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r = (int)(255);
						reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g = (int)(255);
						reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b = (int)(255);

						diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r = (int)(255);
						diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g = (int)(255);
						diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b = (int)(255);

						emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r = (int)(255);
						emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g = (int)(255);
						emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b = (int)(255);
					}


				}
			}
			//gather Anti Aliasing
			colorInt ap = { 0, 0, 0 }, re = { 0, 0, 0 }, diff = { 0, 0, 0 }, emis = { 0, 0, 0 };
			for (int iAA = 0; iAA < AA; iAA++){
				for (int jAA = 0; jAA < AA; jAA++){
					ap.r += apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r;
					ap.g += apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g;
					ap.b += apAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b;

					re.r += reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r;
					re.g += reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g;
					re.b += reAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b;

					diff.r += diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r;
					diff.g += diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g;
					diff.b += diffAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b;

					emis.r += emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].r;
					emis.g += emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].g;
					emis.b += emisAA[omp_get_thread_num()*AA*AA+iAA*AA + jAA].b;
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

	saveData(tempData, width, height);

//	std::cout << "countBIEN, countMAL, total " << countBIEN << " " << countMAL << " " << width*height << "\n";
	std::cout << "Thermography generated successfully \n";
	std::cout << "You can close this window \n";

	//SAVE IMAGE
	std::stringstream ss1;
	ss1 << "../results/real.png";

	std::stringstream ss2;
	ss2 << "../results/apparent.png";

	std::stringstream ss3;
	ss3 << "../results/refl.png";

	std::stringstream ss4;
	ss4 << "../results/emis.png";

	saveImage(reColors, width, height, ss1.str());
	saveImage(apColors, width, height, ss2.str());
	saveImage(diffColors, width, height, ss3.str());
	saveImage(emisColors, width, height, ss4.str());

	delete reColors;
	delete apColors;
	delete diffColors;
};

float* loadSkyTemp(string file){
	float* tsky = (float*)malloc(sizeof(float)*10);

	std::ifstream ifs(file);

	for (int j = 0; j < 10; j++){
		ifs >> tsky[j];
	}
	std::cout << "Sky temps loaded successfully...\n";

	return tsky;
}

int main(){
	std::vector<glm::vec3> sc_vertices;
	std::vector<int> sc_triangles;
	std::vector<int> sc_quads;
	std::vector<int> matIDs;
	std::vector<float> temps;
	settings s = loadSettings("..\\viewSettings");

	load_UCD(("..\\" + s.sceneFile), sc_vertices, sc_triangles, sc_quads, matIDs, temps);
	buildSceneEmbree(sc_vertices, sc_triangles, sc_quads, matIDs, temps);

	loadColormapFromFile(("..\\" + s.colormapFile));
	float* tsky = loadSkyTemp(("..\\" + s.skyTempsFile));
	tmin = s.tmin; tmax = s.tmax;
	tmin_reflected = s.tmin_reflected; tmax_reflected = s.tmax_reflected;

	NRAYS_GLOSSY = s.reflSamples;

	material* matProps = loadMaterials("..\\materials");
	//printMaterials(matProps);

	generateThermography(tsky, matIDs, s,matProps);
	int i;
	cin >> i;
	return 0;
}