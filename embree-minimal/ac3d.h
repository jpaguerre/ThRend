#ifndef AC3D_LOADER
#define AC3D_LOADER

#include <string>
#include <vector>
#include <glm.hpp>

namespace ac3d {

	int load_scene(std::string filename, std::vector<glm::vec3> &sc_vertices, std::vector<int> &sc_triangles, std::vector<int> &sc_quads);
			

}

#endif

