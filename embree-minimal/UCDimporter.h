#ifndef UCDimporter
#define UCDimporter

#include <string>
#include <vector>
#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
using namespace glm;


int load_UCD(std::string filename, 
std::vector<vec3> &sc_vertices, 
std::vector<int> &sc_triangles, 
std::vector<int> &sc_quads, 
std::vector<int> &matIDs,
std::vector<float> &temps){
	ifstream file(filename.c_str());
	cout << "Loading file " << filename << " \n";

	if (!file) 
		return -1;

	string line;
	getline(file, line);

	int nPoints; int nFaces;
	stringstream   linestream(line);
	linestream >> nPoints >> nFaces;
	cout << "Total points: " << nPoints << " \n";
	cout << "Total faces: " <<  nFaces << " \n";
	//load vertices
	for (int i = 0; i < nPoints; i++){
		getline(file, line);
		int id; float x, y, z;
		stringstream   linestream(line);
		linestream >> id >> x >> y >> z;
		sc_vertices.push_back(vec3(x, y, z));
	}
	//load faces
	for (int i = 0; i < nFaces; i++){
		getline(file, line);
		int id,mat; 
		string poly;
		int idx1, idx2, idx3, idx4;
		stringstream   linestream(line);
		linestream >> id >> mat >> poly;
		matIDs.push_back(mat);

		if (poly == "tri"){
			linestream >> idx1 >> idx2 >> idx3;
			sc_triangles.push_back(idx1-1);
			sc_triangles.push_back(idx2-1);
			sc_triangles.push_back(idx3-1);
		}
		else if (poly == "quad"){
			linestream >> idx1 >> idx2 >> idx3 >> idx4;
			sc_quads.push_back(idx1-1);
			sc_quads.push_back(idx2-1);
			sc_quads.push_back(idx3-1);
			sc_quads.push_back(idx4-1);
		}
	}
	cout << "Geometry loaded successfully \n";
	//discard two lines
	getline(file, line);
	getline(file, line);
	//load temperatures
	for (int i = 0; i < nPoints; i++){
		getline(file, line);
		int id; float t;
		stringstream   linestream(line);
		linestream >> id >> t;
		temps.push_back(t);
	}
	cout << "Nodal temperatures loaded successfully \n";

	return 1;
}


#endif

