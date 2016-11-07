#include "grid.h"

double Lattice::Interpolate(const Vector3i &cell, const Vector3d &coords) const
{
	int k1 = cell.z * size_xy + cell.y * size.x + cell.x; // (x, y, z)
	int k2 = k1 + size_xy;                                // (x, y, z + 1)
	int k3 = k1 + size.x;                                 // (x, y + 1, z)
	int k4 = k3 + size_xy;                                // (x, y + 1, z + 1)

	double xd = coords.x - cell.x;
	double yd = coords.y - cell.y;
	double zd = coords.z - cell.z;

	double xd_inv = 1.0 - xd;
	double yd_inv = 1.0 - yd;

	double c00 = data[k1] * xd_inv + data[k1 + 1] * xd;
	double c01 = data[k2] * xd_inv + data[k2 + 1] * xd;
	double c10 = data[k3] * xd_inv + data[k3 + 1] * xd;
	double c11 = data[k4] * xd_inv + data[k4 + 1] * xd;

	double c0 = c00 * yd_inv + c10 * yd;
	double c1 = c01 * yd_inv + c11 * yd;

	return c0 * (1.0 - zd) + c1 * zd;
}

YeeGrid::YeeGrid(const Vector3d &vmin, const Vector3d &vmax, const Vector3i &numCells) :
	vmin(vmin),
	vmax(vmax),
	numCells(numCells),
	cellSize((vmax - vmin) / numCells),

	Ex(numCells + Vector3i(2, 1, 1)),
	Ey(numCells + Vector3i(1, 2, 1)),
	Ez(numCells + Vector3i(1, 1, 2)),

	Bx(numCells + Vector3i(1, 2, 2)),
	By(numCells + Vector3i(2, 1, 2)),
	Bz(numCells + Vector3i(2, 2, 1))
{ }

FieldPoint YeeGrid::Interpolate(const Vector3d &coords) const
{
	Vector3d pos = (coords - vmin) / cellSize;
	Vector3i cell = floor(pos);
	Vector3i cell2 = floor(pos + Vector3d(0.5));

	if (cell.x == numCells.x) cell.x -= 1;
	if (cell.y == numCells.y) cell.y -= 1;
	if (cell.z == numCells.z) cell.z -= 1;
	if (cell2.x == numCells.x) cell2.x -= 1;
	if (cell2.y == numCells.y) cell2.y -= 1;
	if (cell2.z == numCells.z) cell2.z -= 1;

	Vector3i cell_Ex(cell), cell_Ey(cell), cell_Ez(cell);
	cell_Ex.x = cell2.x;
	cell_Ey.y = cell2.y;
	cell_Ez.z = cell2.z;

	Vector3i cell_Bx(cell2), cell_By(cell2), cell_Bz(cell2);
	cell_Bx.x = cell.x;
	cell_By.y = cell.y;
	cell_Bz.z = cell.z;

	FieldPoint f;

	f.E = Vector3d(
		Ex.Interpolate(cell_Ex, pos + Vector3d(0.5, 0, 0)),
		Ey.Interpolate(cell_Ey, pos + Vector3d(0, 0.5, 0)),
		Ez.Interpolate(cell_Ez, pos + Vector3d(0, 0, 0.5))
	);

	f.B = Vector3d(
		Bx.Interpolate(cell_Bx, pos + Vector3d(0, 0.5, 0.5)),
		By.Interpolate(cell_By, pos + Vector3d(0.5, 0, 0.5)),
		Bz.Interpolate(cell_Bz, pos + Vector3d(0.5, 0.5, 0))
	);

	return f;
}

void WeightCurrents(const std::vector<Particle> &particles)
{

}