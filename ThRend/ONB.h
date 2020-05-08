//THIS FILE WAS WRITTEN BY IGNACIO DECIA

#ifndef OrtoNormalBase
#define OrtoNormalBase

#include <glm.hpp>


class ONB {

public:
	ONB(glm::vec3 normal);
	glm::vec3 WorldToLocal(const glm::vec3 &v);
	glm::vec3 LocalToWorld(const glm::vec3 &v);

private:
	glm::vec3 n;
	glm::vec3 s;
	glm::vec3 t;
};


ONB::ONB(glm::vec3 normal)	{
     n = normal;
	if (fabs(n.x) > fabs(n.z)) {
		s.x = -n.y;
		s.y = n.x;
		s.z = 0;
	} else {

		s.x = 0;
		s.y = -n.z;
		s.z = n.y;
	}

	s = normalize(s);
	t = cross(n, s);
}


glm::vec3 ONB::WorldToLocal(const glm::vec3 &v) {
	return glm::vec3(dot(v, s), dot(v, t), dot(v, n));
}


glm::vec3 ONB::LocalToWorld(const glm::vec3 &v) {
	return  (v.x * s) + (v.y * t) + (v.z * n);
}

#endif