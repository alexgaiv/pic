#ifndef _UTILS_H_
#define _UTILS_H_

#define __CL_ENABLE_EXCEPTIONS

#include <CL\cl.hpp>
#include <fstream>

using namespace cl;
using namespace std;

inline Program CreateProgramFromFile(const char *filename, const Context &ctx)
{
	ifstream file(filename);
	if (!file) throw Error(0, "LoadProgramFromFile");

	vector<string> lines;
	string line;
	while (getline(file, line)) {
		lines.push_back(line);
	}

	Program::Sources source;
	for (int i = 0, n = lines.size(); i < n; i++) {
		source.push_back(make_pair(lines[i].c_str(), lines[i].size()));
	}

	return Program(ctx, source);
}

#endif // _UTILS_H_