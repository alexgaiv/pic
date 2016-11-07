#include "grid.h"

double Lattice::Interpolate(const Vector3i &cell, const Vector3d &coords) const
{
	int k1 = cell.z * size_xy + cell.y * size.x + cell.x; // (x, y, z)
	int k2 = k1 + size_xy;                                // (x, y, z + 1)
	int k3 = k1 + size.x;                                 // (x, y + 1, z)
	int k4 = k3 + size_xy;                                // (x, y + 1, z + 1)

	Vector3d c = coords - cell;
	Vector3d c_inv = Vector3d(1.0) - c;

	double c00 = data[k1] * c_inv.x + data[k1 + 1] * c.x;
	double c01 = data[k2] * c_inv.x + data[k2 + 1] * c.x;
	double c10 = data[k3] * c_inv.x + data[k3 + 1] * c.x;
	double c11 = data[k4] * c_inv.x + data[k4 + 1] * c.x;

	double c0 = c00 * c_inv.y + c10 * c.y;
	double c1 = c01 * c_inv.y + c11 * c.y;

	return c0 * c_inv.z + c1 * c.z;
}

void Lattice::Deposit(const Vector3i &cell, const Vector3d &coords, double value)
{
	int k1 = cell.z * size_xy + cell.y * size.x + cell.x; // (x, y, z)
	int k2 = k1 + size_xy;                                // (x, y, z + 1)
	int k3 = k1 + size.x;                                 // (x, y + 1, z)
	int k4 = k3 + size_xy;                                // (x, y + 1, z + 1)

	Vector3d c = coords - cell;
	Vector3d c_inv = Vector3d(1.0) - c;

	data[k1]     += c_inv.x * c_inv.y * c_inv.z * value;
	data[k2]     += c_inv.x * c_inv.y * c.z * value;
	data[k3]     += c_inv.x * c.y * c_inv.z * value;
	data[k4]     += c_inv.x * c.y * c.z * value;
	data[k1 + 1] += c.x * c_inv.y * c_inv.z * value;
	data[k2 + 1] += c.x * c_inv.y * c.z * value;
	data[k3 + 1] += c.x * c.y * c_inv.z * value;
	data[k4 + 1] += c.x * c.y * c.z * value;
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
	Bz(numCells + Vector3i(2, 2, 1)),

	Jx(numCells + Vector3i(2, 1, 1)),
	Jy(numCells + Vector3i(1, 2, 1)),
	Jz(numCells + Vector3i(1, 1, 2)),

	shift_JEx(0.5, 0, 0),
	shift_JEy(0, 0.5, 0),
	shift_JEz(0, 0, 0.5),

	shift_Bx(0, 0.5, 0.5),
	shift_By(0.5, 0, 0.5),
	shift_Bz(0.5, 0.5, 0)
{ }

FieldPoint YeeGrid::InterpolateField(const Vector3d &coords) const
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
		Ex.Interpolate(cell_Ex, pos + shift_JEx),
		Ey.Interpolate(cell_Ey, pos + shift_JEy),
		Ez.Interpolate(cell_Ez, pos + shift_JEz)
	);

	f.B = Vector3d(
		Bx.Interpolate(cell_Bx, pos + shift_Bx),
		By.Interpolate(cell_By, pos + shift_By),
		Bz.Interpolate(cell_Bz, pos + shift_Bz)
	);

	return f;
}

void YeeGrid::DepositCurrents(const Particle &pt)
{
	Vector3d j = pt.factor * pt.charge * pt.Velocity();

	Vector3d pos = (pt.coords - vmin) / cellSize;
	Vector3i cell = floor(pos);
	Vector3i cell2 = floor(pos + Vector3d(0.5));

	if (cell.x == numCells.x) cell.x -= 1;
	if (cell.y == numCells.y) cell.y -= 1;
	if (cell.z == numCells.z) cell.z -= 1;
	if (cell2.x == numCells.x) cell2.x -= 1;
	if (cell2.y == numCells.y) cell2.y -= 1;
	if (cell2.z == numCells.z) cell2.z -= 1;

	Vector3i cell_Jx(cell), cell_Jy(cell), cell_Jz(cell);
	cell_Jx.x = cell2.x;
	cell_Jy.y = cell2.y;
	cell_Jz.z = cell2.z;

	Jx.Deposit(cell_Jx, pos + shift_JEx, j.x);
	Jy.Deposit(cell_Jy, pos + shift_JEy, j.y);
	Jz.Deposit(cell_Jz, pos + shift_JEz, j.z);
}