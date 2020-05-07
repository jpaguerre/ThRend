#ifndef MATERIALimporter
#define MATERIALimporter

#include <string>
#include <vector>
#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>


using namespace std;
using namespace glm;

typedef struct{
	string name;
	int UCD_id;
	float normal_emissivity;
	float diffuse_fraction;
	float specular_lobe_size;
	float * emisTable;
	bool custom;
} material;

vector<int> ids;

double deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

float * computeEmissivityCurve(material m){
	float* emisTable = (float*)malloc(sizeof(float) * 91);
	float e = m.normal_emissivity;
	float df = m.diffuse_fraction;
	float r = 1.0 - e;
	for (int i = 0; i < 91; i++){
		float teta = deg2rad(90 - i);
		float schlick = 1 - (r + (1 - r)*(pow(1 - cos(teta), 5.0)));
		emisTable[i] = (df*e) + (1 - df)*schlick;
	}
	return emisTable;
}


material* loadMaterials(std::string filename){
	material* matProps = (material*)malloc(sizeof(material) * 64);
	material m; m.UCD_id = -1; m.custom = false;
	ifstream file(filename.c_str());

	cout << "Loading file " << filename << " \n";
	string line;
	while (getline(file, line)){
		stringstream   linestream(line);
		string id;
		linestream >> id;
		if (id.size() < 0 || (id.size() > 0 && id.at(0) == '#')){
			
		}
		else if (id == "name"){
			//allocate new material
			if (m.UCD_id != -1){
				ids.push_back(m.UCD_id);
				if (!m.custom)
					m.emisTable = computeEmissivityCurve(m);
				memcpy(matProps + (m.UCD_id), &m, sizeof(material));
				m.custom = false;
			}
			string x;
			linestream >> x;
			m.name = x;
		}
		else if (id == "UCD_id"){
			int x;
			linestream >> x;
			m.UCD_id = x;
		}
		else if (id == "normal_emissivity"){
			float x;
			linestream >> x;
			m.normal_emissivity = x;
		}
		else if (id == "diffuse_fraction"){
			float x;
			linestream >> x;
			m.diffuse_fraction = x;
		}
		else if (id == "specular_lobe_size"){
			float x;
			linestream >> x;
			m.specular_lobe_size = x;
		}
		else if (id == "emissivity_curve"){
			float x;
			m.custom = true;
			m.emisTable = (float*)malloc(sizeof(float) * 91);
			for (int i = 0; i < 91; i++){
				linestream >> x;
				m.emisTable[i] = x;
			}
		}
	}
	//copy last one
	if (m.UCD_id != -1){
		ids.push_back(m.UCD_id);
		if (!m.custom)
			m.emisTable = computeEmissivityCurve(m);
		memcpy(matProps + (m.UCD_id), &m, sizeof(material));
		m.custom = false;
	}
	cout << "Materials loaded succesfully \n";
	return matProps;
}

void printMaterials(material* matProps){
	for (int i = 0; i < ids.size(); i++){
		material m = matProps[ids[i]];
		cout << "name " << m.name << "\n";
		cout << "UCD_id " << m.UCD_id << "\n";
		if (!m.custom){
			cout << "normal_emissivity " << m.normal_emissivity << "\n";
			cout << "diffuse_fraction " << m.diffuse_fraction << "\n";
		}
		cout << "specular_lobe_size " << m.specular_lobe_size << "\n";

		for (int j = 0; j < 91; j++)
			cout << m.emisTable[j] << " ";
		cout << "\n";
	}
}


#endif

