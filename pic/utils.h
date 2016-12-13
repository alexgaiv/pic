#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>
#include <math.h>

#ifndef PRESICION
#define PRESICION 0.0001;
#endif

#define FOR3(i, j, k, s) \
	for (int i = 0; i < s.x; i++) \
	for (int j = 0; j < s.y; j++) \
	for (int k = 0; k < s.z; k++)

inline double frand(double min, double max)
{
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

inline bool CmpReal(double a, double b) {
	return fabs(a - b) <= PRESICION;
}

#endif // _UTILS_H_