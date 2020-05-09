#ifndef SETimporter
#define SETimporter

#include <string>
#include <vector>
#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
using namespace glm;

typedef struct{
	string sceneFile;
	string skyTempsFile;
	string colormapFile;
	vec3 cameraCenter;
	vec3 cameraDirection;
	vec3 cameraUp;
	float fovVertical;
	int	imageWidth;
	int	imageHeight;
	int aa;
	int MAX_BOUNCES;
	int reflSamples;
	float tmin;
	float tmax;
	float tmin_reflected;
	float tmax_reflected;
} settings;


settings loadSettings(std::string filename){
	ifstream file(filename.c_str());
	cout << "Loading file " << filename << " \n";
	string line;
	settings s;
	while (getline(file, line)){
		stringstream   linestream(line);
		string id;
		linestream >> id;
		if (id.size() < 0 || (id.size() > 0 && id.at(0) == '#')){
			
		}
		else if (id == "sceneFile"){
			string x;
			linestream >> x;
			s.sceneFile = x;
		}
		else if (id == "skyTempsFile"){
			string x;
			linestream >> x;
			s.skyTempsFile = x;
		}
		else if (id == "colormapFile"){
			string x;
			linestream >> x;
			s.colormapFile = x;
		}
		else if (id == "cameraCenter"){
			float x, y, z;
			linestream >> x >> y >> z;
			s.cameraCenter = vec3(x, y, z);
		}
		else if (id == "cameraDirection"){
			float x, y, z;
			linestream >> x >> y >> z;
			s.cameraDirection = vec3(x, y, z);
		}
		else if (id == "cameraUp"){
			float x, y, z;
			linestream >> x >> y >> z;
			s.cameraUp = vec3(x, y, z);
		}
		else if (id == "fovVertical"){
			float x;
			linestream >> x ;
			s.fovVertical = x;
		}
		else if (id == "imageWidth"){
			int x;
			linestream >> x;
			s.imageWidth = x;
		}
		else if (id == "imageHeight"){
			int x;
			linestream >> x;
			s.imageHeight = x;
		}
		else if (id == "aa"){
			int x;
			linestream >> x;
			s.aa = x;
		}
		else if (id == "MAX_BOUNCES"){
			int x;
			linestream >> x;
			s.MAX_BOUNCES = x;
			if (s.MAX_BOUNCES < 1){
				std::cout << "ERROR: please input MAX_BOUNCES greater than 0.\n";
				s.MAX_BOUNCES = 1;
			}
		}
		else if (id == "reflSamples"){
			int x;
			linestream >> x;
			s.reflSamples = x;
		}
		else if (id == "tmin"){
			float x;
			linestream >> x;
			s.tmin = x;
		}
		else if (id == "tmax"){
			float x;
			linestream >> x;
			s.tmax = x;
		}
		else if (id == "tmin_reflected"){
			float x;
			linestream >> x;
			s.tmin_reflected = x;
		}
		else if (id == "tmax_reflected"){
			float x;
			linestream >> x;
			s.tmax_reflected = x;
		}
	}
	cout << "View settings loaded succesfully \n";
	return s;
}

void printSettings(settings s){
	cout << "cameraCenter " << s.cameraCenter.x << " " << s.cameraCenter.y << " " << s.cameraCenter.z << "\n";
	cout << "cameraDirection " << s.cameraDirection.x << " " << s.cameraDirection.y << " " << s.cameraDirection.z << "\n";
	cout << "cameraUp " << s.cameraUp.x << " " << s.cameraUp.y << " " << s.cameraUp.z << "\n";
	cout << "fovVertical " << s.fovVertical << "\n";
	cout << "imageWidth " << s.imageWidth << "\n";
	cout << "imageHeight " << s.imageHeight << "\n";
	cout << "aa " << s.aa << "\n";
}


#endif

