#define _USE_MATH_DEFINES
#define GLM_ENABLE_EXPERIMENTAL
#define TASKING_TBB
//#define EMBREE_BACKFACE_CULLING 1


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

#include "emissivity.h"
#include "UCDimporter.h"
#include "colormap.h"
#include "reflections.h"
#include "ONB.h"
#include <chrono>  // for high_resolution_clock

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
 
#include <windows.h>
#include <process.h> 
#include "FreeImage.h"
#include "EmbreeUtils.h"
#include "settingsLoader.h"
#include "materialLoader.h"

#include <fstream>
#include <iostream>

const float STE_BOLZ = 5.670367e-08;
const int NUMTHREADS = 8;


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

void printProgress(double percentage) {
	if (percentage > 0.95) percentage=1.0;
	int barLength = 60;	float pos = percentage * barLength;
	std::cout << " [";
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

	//allocate memory for the output images
	RGBQUAD *apColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *reColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *diffColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	RGBQUAD *emisColors = (RGBQUAD*)malloc(width*height*sizeof(RGBQUAD));
	float *tempData = (float*)malloc(width*height*sizeof(float));

	int countMAL = 0;
	int countBIEN = 0;
	totalProcessed = 0;

	//Beckers subdivision for the case of diffuse reflection
	float rings = ceil((-1 + sqrt(NRAYS_GLOSSY)) / 2);
	int nRaysDiffuse = 1 + 4 * rings*(rings + 1); // n >= al n original
	float	Drho = 1 / (2 * rings + 1); //Delta rho
	//
	omp_set_num_threads(NUMTHREADS);
	float* reflTAnt = (float*)malloc(NUMTHREADS*sizeof(float));
	//allocate memory for the intermediate antialisiasing results
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
							if (specN == -1)//DIFFUSE REFLECTION
								reflT = getDiffuselyReflectedTemperature(tsky, orig, base, rings, Drho, nRaysDiffuse);
							else{           //GLOSSY REFLECTION
								reflT = getGlossyReflectedTemperature(tsky, orig, originalDir, hitNormal, angulo, specN, gID, pID, reflTAnt);
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

	std::cout << "Thermography generated and saved successfully \n";
	std::cout << "You can close this window \n";

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