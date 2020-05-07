#include "ac3d.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

const float emission_factor = 5;

glm::vec3 read_vertex(ifstream& ifs) {
	glm::vec3 v;
	string line;
	getline(ifs, line);
	stringstream ostr(line);
	ostr >> v.x >> v.y >> v.z;

	if (abs(v.z - -0.3)< 1e-3 ){
		float delta=1*((v.y - 80)/-127.);

		v.z = 0.01 + delta;
	}
	return v;
}

std::vector<int> read_surface(ifstream& ifs, unsigned int vert_offset) {
	std::vector<int> s;

	while(1) {
		string line;
		getline(ifs, line);
		stringstream ostr(line);
		string token;
		ostr >> token;

		if(token == "refs") {
			int refs;
			ostr >> refs;
			for(int i=0; i<refs; ++i) {
				string line;
				getline(ifs, line);
				stringstream ostr(line);
				int v;
				ostr >> v;
				s.push_back(v+vert_offset);
			}
			break;
		}
	}
	return s;
}

void read_object(ifstream& ifs, std::vector<glm::vec3> &sc_vertices, std::vector<int> &sc_triangles, std::vector<int> &sc_quads) {
	/*
	OBJECT %s
	*name %s
	*data %d
	*texture %s
	*texrep %f %f
	*rot %f %f %f  %f %f %f  %f %f %f
	*loc %f %f %f
	*url %s
	*numvert %d
	numvert lines of %f %f %f
	*numsurf %d
	*SURF %d
	*mat %d
	refs %d
	refs lines of %d %f %f
	kids %d 
	*/
	size_t vert_offset = sc_vertices.size();
    
//    nodehierarchy nh;
	unsigned int indexnodes;
    while(1) {
		string line;
		getline(ifs, line);
		stringstream ostr(line);
		string token;
		ostr >> token;

		if(token == "numvert") {
			int numvert;
			ostr >> numvert;
			for(int i=0; i<numvert; ++i)
				sc_vertices.push_back(read_vertex(ifs));
		} else if(token == "numsurf") {
			unsigned int numsurf;
			ostr >> numsurf;
			for(unsigned int i=0; i<numsurf; ++i) {
				std::vector<int> surf;
                surf=read_surface(ifs, static_cast<unsigned int>(vert_offset));
				if (surf.size()==3)
					for (int j = 0; j < surf.size();j++)
						sc_triangles.push_back(surf[j]);
				if (surf.size() == 4)
					for (int j = 0; j < surf.size(); j++)
						sc_quads.push_back(surf[j]);
				/*if (surf.size() == 4){
					sc_triangles.push_back(surf[0]);
					sc_triangles.push_back(surf[1]);
					sc_triangles.push_back(surf[2]);

					sc_triangles.push_back(surf[0]);
					sc_triangles.push_back(surf[2]);
					sc_triangles.push_back(surf[3]);				
				}*/				
			}
		} else if(token == "kids")
			break;
	}
}

int ac3d::load_scene(std::string filename, std::vector<glm::vec3> &sc_vertices, std::vector<int> &sc_triangles, std::vector<int> &sc_quads)
{
    ifstream ifs(filename.c_str());

	if(!ifs)
		return -1;

	string line;
	getline(ifs, line);

	if(line != "AC3Db")
		return -1;

	while(!ifs.eof()) {
		getline(ifs, line);
		stringstream ostr(line);
		string token;
		ostr >> token;

		if(token == "OBJECT") {
			ostr >> token;
			if(token == "world")
				getline(ifs, line);
			else if(token == "poly")
				read_object(ifs, sc_vertices, sc_triangles, sc_quads);
		}
	}
	
	return 0;
}


