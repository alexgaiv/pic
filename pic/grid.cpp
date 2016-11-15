#define _USE_MATH_DEFINES
#include "grid.h"
#include <math.h>

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

YeeGrid::YeeGrid(const Vector3d &vmin, const Vector3d &vmax, const Vector3i &numInnerCells) :
	vmin(vmin),
	vmax(vmax),
	numCells(numInnerCells + Vector3i(2)),
	cellSize((vmax - vmin) / numInnerCells),

	Ex(numCells + Vector3i(0, 1, 1)),
	Ey(numCells + Vector3i(1, 0, 1)),
	Ez(numCells + Vector3i(1, 1, 0)),

	Bx(numCells + Vector3i(1, 0, 0)),
	By(numCells + Vector3i(0, 1, 0)),
	Bz(numCells + Vector3i(0, 0, 1)),

	Jx(numCells + Vector3i(0, 1, 1)),
	Jy(numCells + Vector3i(1, 0, 1)),
	Jz(numCells + Vector3i(1, 1, 0)),

	shift_JEx(0.5, 1, 1),
	shift_JEy(1, 0.5, 1),
	shift_JEz(1, 1, 0.5),

	shift_Bx(1, 0.5, 0.5),
	shift_By(0.5, 1, 0.5),
	shift_Bz(0.5, 0.5, 1)
{ }

FieldPoint YeeGrid::InterpolateField(const Vector3d &coords) const
{
	Vector3d pos = (coords - vmin) / cellSize;
	Vector3i cell = floor(pos);
	Vector3i cell2 = floor(pos + Vector3d(0.5));

	Vector3i cell_Ex(cell + Vector3i(1));
	Vector3i cell_Ey(cell_Ex);
	Vector3i cell_Ez(cell_Ex);
	cell_Ex.x = cell2.x;
	cell_Ey.y = cell2.y;
	cell_Ez.z = cell2.z;

	Vector3i cell_Bx(cell2);
	Vector3i cell_By(cell2);
	Vector3i cell_Bz(cell2);
	cell_Bx.x = cell.x + 1;
	cell_By.y = cell.y + 1;
	cell_Bz.z = cell.z + 1;

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

	Vector3i cell_Jx(cell + Vector3i(1));
	Vector3i cell_Jy(cell_Jx);
	Vector3i cell_Jz(cell_Jx);
	cell_Jx.x = cell2.x;
	cell_Jy.y = cell2.y;
	cell_Jz.z = cell2.z;

	Jx.Deposit(cell_Jx, pos + shift_JEx, j.x);
	Jy.Deposit(cell_Jy, pos + shift_JEy, j.y);
	Jz.Deposit(cell_Jz, pos + shift_JEz, j.z);
}

void YeeGrid::SolveField(double dt)
{
	double cdt = c * dt;
	Vector3d d = Vector3d(1.0) / cellSize;
	double k = 4 * M_PI * dt;

	// Ex
	for (int i = 0; i < Ex.GetSize().x; i++)
		for (int j = 0; j < Ex.GetSize().y; j++)
			for (int k = 0; k < Ex.GetSize().z; k++)
			{
				int j1 = j;
				int j2 = j - 1;
				int k1 = k;
				int k2 = k - 1;

				if (j == 0)
					j2 = Bz.GetSize().y - 1;
				else if (j == Ex.GetSize().y - 1)
					j1 = 0;

				if (k == 0)
					k2 = By.GetSize().z - 1;
				else if (k == Ex.GetSize().z - 1)
					k1 = 0;

				Ex(i, j, k) += cdt *
					((Bz(i, j1, k) - Bz(i, j2, k)) * d.y -
					(By(i, j, k1) - By(i, j, k2)) * d.z) - 
					k * Jx(i, j, k);
			}

	// Ey
	for (int i = 0; i < Ey.GetSize().x; i++)
		for (int j = 0; j < Ey.GetSize().y; j++)
			for (int k = 0; k < Ey.GetSize().z; k++)
			{
				int i1 = i;
				int i2 = i - 1;
				int k1 = k;
				int k2 = k - 1;

				if (i == 0)
					i2 = Bz.GetSize().x - 1;
				else if (i == Ey.GetSize().x - 1)
					i1 = 0;

				if (k == 0)
					k2 = Bx.GetSize().z - 1;
				else if (k == Ey.GetSize().z - 1)
					k1 = 0;

				Ey(i, j, k) -= cdt *
					((Bz(i1, j, k) - Bz(i2, j, k)) * d.x -
					(Bx(i, j, k1) - Bx(i, j, k2)) * d.z) -
					k * Jy(i, j, k);
			}

	// Ez
	for (int i = 0; i < Ez.GetSize().x; i++)
		for (int j = 0; j < Ez.GetSize().y; j++)
			for (int k = 0; k < Ez.GetSize().z; k++)
			{
				int i1 = i;
				int i2 = i - 1;
				int j1 = j;
				int j2 = j - 1;

				if (i == 0)
					i2 = By.GetSize().x - 1;
				else if (i == Ez.GetSize().x - 1)
					i1 = 0;

				if (j == 0)
					j2 = Bx.GetSize().y - 1;
				else if (k == Ez.GetSize().z - 1)
					j1 = 0;

				Ez(i, j, k) += cdt *
					((By(i1, j, k) - By(i2, j, k)) * d.x -
					(Bx(i, j1, k) - Bx(i, j2, k)) * d.y) -
					k * Jz(i, j, k);
			}

	// Bx
	for (int i = 0; i < Bx.GetSize().x; i++)
		for (int j = 0; j < Bx.GetSize().y; j++)
			for (int k = 0; k < Bx.GetSize().z; k++)
			{
				Bx(i, j, k) -= cdt *
					((Ez(i, j + 1, k) - Ez(i, j, k)) * d.y -
					(Ey(i, j, k + 1) - Ey(i, j, k)) * d.z);
			}

	// By
	for (int i = 0; i < By.GetSize().x; i++)
		for (int j = 0; j < By.GetSize().y; j++)
			for (int k = 0; k < By.GetSize().z; k++)
			{
				By(i, j, k) += cdt *
					((Ez(i + 1, j, k) - Ez(i, j, k)) * d.x -
					(Ex(i, j, k + 1) - Ex(i, j, k)) * d.z);
			}

	// Bz
	for (int i = 0; i < Bz.GetSize().x; i++)
		for (int j = 0; j < Bz.GetSize().y; j++)
			for (int k = 0; k < Bz.GetSize().z; k++)
			{
				Bz(i, j, k) -= cdt *
					((Ey(i + 1, j, k) - Ey(i, j, k)) * d.x -
					(Ex(i, j + 1, k) - Ex(i, j, k)) * d.y);
			}
}