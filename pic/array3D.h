#ifndef _ARRAY_3D_H_
#define _ARRAY_3D_H_

#include <vector>
#include "mathtypes.h"

template<class T>
class Array3D
{
public:
	Array3D(const Vector3i &size) :
		size(size),
		data(size.x * size.y * size.z)
	{ }

	Vector3i GetSize() const { return size; }

	T &operator()(int i, int j, int k) {
		return data[(k * size.y + j) * size.x + i];
	}
	const T &operator()(int i, int j, int k) const {
		return data[(k * size.y + j) * size.x + i];
	}
private:
	Vector3i size;
	std::vector<T> data;
};

#endif // _ARRAY_3D_H_