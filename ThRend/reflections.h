#ifndef REFLECTIONS
#define REFLECTIONS

#include <string>
#include <vector>
#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include "EmbreeUtils.h"
#include "materialLoader.h"
#include "ONB.h"
#include <omp.h>

using namespace std;
using namespace glm;

int NRAYS_GLOSSY;
int MAX_BOUNCES;
inline float power4(float t){ return t*t*t*t; }

Result followSpecularPath(float* tsky, glm::vec3 &reflDir, glm::vec3 &hitPoint, std::vector<int> &matIDs, material *matProps){
	float rayCon = 1;
	float flux = 0;
	Result r;
	float t;
	float emis;
	int globalID;
	int bounces = 0;
	while ((rayCon > 0.001) && (bounces<MAX_BOUNCES)){
		bounces++;
		RTCRay ray;
		ray.org_x = hitPoint.x; ray.org_y = hitPoint.y;	ray.org_z = hitPoint.z;
		ray.tfar = 1000.0f; ray.tnear = EPS; ray.time = 0; ray.id = 0;
		//avoid self hit
		ray.flags = -1;
		ray.dir_x = reflDir.x; 	ray.dir_y = reflDir.y; ray.dir_z = reflDir.z;

		RTCRayHit query;
		RTCIntersectContext context;
		rtcInitIntersectContext(&context);
		query.ray = ray;
		query.hit.geomID = RTC_INVALID_GEOMETRY_ID;
		query.hit.primID = RTC_INVALID_GEOMETRY_ID;
		rtcIntersect1(eScene, &context, &query);
		
		if (query.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
			glm::vec3 color(0.0f, 0.0f, 0.0f);

			r.gID = query.hit.geomID;
			r.pID = query.hit.primID;
			float u = query.hit.u;
			float v = query.hit.v;

			globalID = r.pID;
			if (r.gID == quadID)
				globalID += triOffset;

			{
				rtcInterpolate0(rtcGetGeometry(eScene, r.gID), r.pID, u, v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &t, 1);
			}
			// getSpecularDirection
			glm::vec3 originalDir = glm::vec3(query.ray.dir_x, query.ray.dir_y, query.ray.dir_z);
			glm::vec3 hitNormal = glm::normalize(glm::vec3(query.hit.Ng_x, query.hit.Ng_y, query.hit.Ng_z));
			reflDir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
			hitPoint = glm::vec3(query.ray.org_x, query.ray.org_y, query.ray.org_z) + query.ray.tfar * originalDir + hitNormal*EPS;
			// getDirectionalEmissivity
			ONB base(hitNormal);
			glm::vec3 dirLocal = base.WorldToLocal(-originalDir);
			dirLocal = glm::normalize(dirLocal);
			float angulo = (asin(dirLocal.z) * 180 / M_PI);
			if (angulo < 0){
				emis = 1.0;
			}
			else{
				int ang1 = floor(angulo);
				int ang2 = ceil(angulo);
				float coef = angulo - ang1;
				float * emiT = matProps[matIDs[globalID]].emisTable;
				emis = (1 - coef)*emiT[ang1] + coef*emiT[ang2];
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
			emis = 1.0;
		}
		flux += rayCon*emis*power4(t);
		rayCon *= (1.0 - emis);
	}
	//si hay un sobrante, se lo asigno todo a la ultima temperatura
	if (rayCon > 0)
		flux += rayCon*power4(t);

	r.val = flux;
	return r;
}


Result followSpecularPath2(float* tsky, glm::vec3 reflDir, glm::vec3 hitPoint, float globalID){
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

float getGlossyReflectedFlux(float refl,float* tsky, vec3 orig, vec3 originalDir, vec3 hitNormal, 
					float angulo, float specN, float* reflTAnt, std::vector<int> &matIDs, material *matProps){
	float reflFlux = 0;

	vec3 dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
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
			vec3 hitPoint = orig;
			Result res = followSpecularPath(tsky, perturbed, hitPoint, matIDs, matProps);
			reflFlux += res.val;
			//Result res = followSpecularPath2(tsky, perturbed, orig, -1);
			//reflFlux += power4(res.val);
			count++;
		}
	if (count > 0){
		reflFlux = reflFlux / (count);
		reflTAnt[omp_get_thread_num()] = reflFlux;
	}
	else{
		reflFlux = reflTAnt[omp_get_thread_num()];
	}
	return reflFlux;
}

float getDiffuselyReflectedTemperature(float* tsky, glm::vec3 orig, ONB base, float rings, float Drho, int nRaysDiffuse){
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
					t = (1 - rest)*tsky[9 - angleL] + rest*tsky[9 - angleU];//- 5.
				}
				else
					t = 283.0;
				reflT += t;
			}
		}
		tilespering = (ring + 1) * 4;
	}
	return reflT / nRaysDiffuse;
}


#endif

