// json.cpp : This file reads json with vectors.
// Author: zz1918@nyu.edu
// JSON parser library (https://github.com/nlohmann/json)
#include <fstream>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "json.hpp"
using json = nlohmann::json;
using namespace Eigen;
using namespace std;

Vector3d read_vec3(const json& x)
{
	return Vector3d(stod(string(x[0])), stod(string(x[1])), stod(string(x[2])));
}

Vector4d read_vec4(const json& x)
{
	return Vector4d(stod(string(x[0])), stod(string(x[1])), stod(string(x[2])), stod(string(x[3])));
}

json read_json(const char* filename)
{
	ifstream in(filename);
	json data;
	in >> data;
	return data;
}