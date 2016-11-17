#ifndef _GRID_H_
#define _GRID_H_

#include "mathtypes.h"
#include "particle.h"
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
	void Deposit(const Vector3i &cell, const Vector3d &coords, double value);
private:
	std::vector<double> data;
	Vector3i size;
	int size_xy;
};

class YeeGrid
{
private:
	Vector3i numCells;
public:
	YeeGrid(const Vector3d &vmin, const Vector3d &vmax, const Vector3i &numInnerCells);

	Vector3d GetMin() const { return vmin; }
	Vector3d GetMax() const { return vmax; }
	Vector3i GetNumCells() const { return numCells; }
	Vector3d GetCellSize() const { return cellSize; }

	FieldPoint InterpolateField(const Vector3d &coords) const;
	void DepositCurrents(const Particle &particle);
	void SolveField(double dt);

	Vector3d
		shift_JEx, shift_JEy, shift_JEz,
		shift_Bx, shift_By, shift_Bz;

	Lattice Ex, Ey, Ez;
	Lattice Bx, By, Bz;
	Lattice Jx, Jy, Jz;
private:
	Vector3d vmin, vmax;
	Vector3d cellSize;

	void pbc_E();
	void pbc_B();
};

#endif // _GRID_H_