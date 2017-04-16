#ifndef _UTILS_H_
#define _UTILS_H_

#include "common.h"
#include "mathtypes.h"
#include <fstream>

#define FOR3(i, j, k, s) \
	for (int i = 0; i < (s).x; i++) \
	for (int j = 0; j < (s).y; j++) \
	for (int k = 0; k < (s).z; k++)

#define FOR3_(i, j, k, s) \
	for (int i = 1; i < (s).x - 1; i++) \
	for (int j = 1; j < (s).y - 1; j++) \
	for (int k = 1; k < (s).z - 1; k++)

inline cl_float3 v2v(const Vector3f &v) {
	cl_float3 vec = { v.x, v.y, v.z };
	return vec;
}

inline cl_double3 v2v(const Vector3d &v) {
	cl_double3 vec = { v.x, v.y, v.z };
	return vec;
}

inline cl_int3 v2v(const Vector3i &v) {
	cl_int3 vec = { v.x, v.y, v.z };
	return vec;
}

inline cl::Program CreateProgramFromFile(const char *filename, const cl::Context &ctx)
{
	std::ifstream file(filename);
	if (!file) throw cl::Error(0, "CreateProgramFromFile");

	std::vector<std::string> lines;
	std::string line;
	while (getline(file, line)) {
		lines.push_back(line + '\n');
	}

	cl::Program::Sources source;
	for (int i = 0, n = (int)lines.size(); i < n; i++) {
		source.push_back(std::make_pair(lines[i].c_str(), lines[i].size()));
	}

	return cl::Program(ctx, source);
}

#endif // _UTILS_H_