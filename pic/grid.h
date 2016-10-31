#ifndef _GRID_H_
#define _GRID_H_

#include "mathtypes.h"
#include <vector>

struct FieldPoint
{
	Vector3d E;
	Vector3d B;
};

class Lattice
{
public:
	Lattice() {  }
	Lattice(const Vector3i &size) :
		data(size.x * size.y * size.z),
		size(size),
		size_xy(size.x * size.y)
	{ }

	Vector3i GetSize() const { return size; }
	std::vector<double> &GetData() { return data; }
	const std::vector<double> &GetData() const { return data; }

	const double &operator()(int i, int j, int k) const {
		return data[(k * size.y + j) * size.x + i];
	}
	double &operator()(int i, int j, int k) {
		return data[(k * size.y + j) * size.x + i];
	}

	double Interpolate(const Vector3i &cell, const Vector3d &coords) const;
private:
	std::vector<double> data;
	Vector3i size;
	int size_xy;
};

class YeeGrid
{
public:
	YeeGrid(const Vector3d &vmin, const Vector3d &vmax, const Vector3i &numCells);

	Vector3d GetMin() const { return vmin; }
	Vector3d GetMax() const { return vmax; }
	Vector3i GetNumCells() const { return numCells; }

	FieldPoint Interpolate(const Vector3d &coords) const;

	Lattice Ex, Ey, Ez;
	Lattice Bx, By, Bz;

private:
	Vector3d vmin, vmax;
	Vector3i numCells;
	Vector3d cellSize;
};

#endif // _GRID_H_