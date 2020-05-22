#ifndef REFLECTIONS
#define REFLECTIONS

#include <string>
#include <vector>
#include <glm.hpp>
#include "glm/gtx/string_cast.hpp"
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

	r.val = power4(t);;
	return r;
}

// Implemented by jpaguerre based on eduardof MATLAB version of:
// [Walter et al] "Microfacet Models for Refraction through Rough Surfaces"
// We use GGX sampling and add weight with G functions
inline float G1vmGGX(vec3 wo, vec3 nn, vec3 mm, float thetaV, float alphaG){
	//eq. 34
	float Xplus = (dot(wo, mm) / dot(wo, nn))>0 ? 1.0 : 0.0;
	float G1om = Xplus * (2.0 / (1 + sqrt(1 + alphaG*alphaG*tan(thetaV)*tan(thetaV))));
	return G1om;
}


vec3 specLobeMicrofacetGGX(ONB baseNorm, vec3 wi, float e1, float e2, float alphaG, float &weight, float rotationAngle){
	//eqs. 35, 36
	vec3 nn(0., 0., 1.);
	float thetam = atan(alphaG*sqrt(e1)/sqrt(1-e1));
	float psim = 2.0f * M_PI*e2;

	float x = sin(thetam)*cos(psim);
	float y = sin(thetam)*sin(psim);
	float z = cos(thetam);
	//rotate around Z axis when AA is ON. If AA is off, rotationAngle=0
	float rotated_x = x*cos(rotationAngle) - y*sin(rotationAngle);
	float rotated_y = x*sin(rotationAngle) + y*cos(rotationAngle);
	//
	vec3 mm(rotated_x, rotated_y, z); mm = normalize(mm);

	wi = baseNorm.WorldToLocal(wi);	wi = normalize(wi);
	vec3 wo = 2.0f * (dot(mm, wi)) * mm - wi; wo = normalize(wo);
	if (wo.z < 0)
		return vec3(0, 0, 0);
	else{
		//eqs. 23,41 
		float G1om = G1vmGGX(wo, nn, mm, acos(dot(wo, nn)), alphaG);
		float G1im = G1vmGGX(wi, nn, mm, acos(dot(wi, nn)), alphaG);
		float Giom = G1im*G1om;
		weight = abs(dot(wi, mm))*Giom / (abs(dot(wi, nn))*abs(dot(mm, nn)));
		wo = baseNorm.LocalToWorld(wo); wo = normalize(wo);
		return wo;
	}
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
				float angulo, float alphaG, float* reflTAnt, std::vector<int> &matIDs, material *matProps, float rotationAngle){
	float reflFlux = 0;

	vec3 dir = 2.0f * glm::dot(hitNormal, -originalDir) * hitNormal + originalDir;
	int nRaysX = round(sqrt(NRAYS_GLOSSY));
	int count = 0;
	float weightTot = 0;
	ONB baseNormal(normalize(hitNormal));
	originalDir = normalize(originalDir);

	for (int rx = 0; rx < nRaysX; rx++)
		for (int ry = 0; ry < nRaysX; ry++){
			vec2 Xi = Hammersley(rx*nRaysX + ry, nRaysX*nRaysX);
			float x = Xi.x;
			float y = Xi.y;
			float weight;
			vec3 perturbed = specLobeMicrofacetGGX(baseNormal, -originalDir, x, y, alphaG, weight, rotationAngle);
			if (length(perturbed) < 0.001){
				//reflected vec is autohit
			}
			else{
				vec3 hitPoint = orig;
				Result res = followSpecularPath(tsky, perturbed, hitPoint, matIDs, matProps);
				reflFlux += res.val * weight;
				weightTot += weight;
				count++;
			}
		}
	if (count > 0){
		reflFlux = reflFlux / (weightTot);
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

