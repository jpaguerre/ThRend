#include <glm.hpp>
#include <math.h>

using namespace glm;

int nRings = 10;
int rings[] = { 1, 6, 12, 18, 23, 28, 32, 35, 37, 38 };
int ringsacum[] = { 1, 7, 19, 37, 60, 88, 120, 155, 192, 230 };
float ringslat[] = { 0.081332, 0.2451, 0.4111, 0.5769, 0.7431, 0.9090, 1.0747, 1.240, 1.405, 1.570 };

static int beckers(const vec3 &v) {

	if (v.z <= 0.0)
		return 0;


	float azim  = atan2(v.y, v.x);
	if (azim < 0)
		azim = 2 * M_PI + azim;
	float zenit = atan2(sqrt(v.x*v.x + v.y*v.y), v.z);
    int l = 0;
	int r = nRings - 1;
	int mid = (int)floorf((l + r) / 2);
	while (l < r) {
		if (zenit > ringslat[mid])
			l = mid + 1;
		else if (zenit < ringslat[mid])
			r = mid;
		else
			break;
		mid = (int)floorf((l + r) / 2);
	}
	int ringid = mid;

	int patch;
	if (ringid == 0)
		patch = 1;
	else {
		float anglestep = 2*M_PI / rings[ringid];
		//cuidado aqui con overflow
		int ringint = (int) ceil(azim / anglestep);
		patch =  ringsacum[ringid-1] + ringint;
	}
	return patch;
}