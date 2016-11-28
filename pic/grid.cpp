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
	//int k1 = cell.z * size_xy + cell.y * size.x + cell.x; // (x, y, z)
	//int k2 = k1 + size_xy;                                // (x, y, z + 1)
	//int k3 = k1 + size.x;                                 // (x, y + 1, z)
	//int k4 = k3 + size_xy;                                // (x, y + 1, z + 1)

	Vector3d c = coords - cell;
	Vector3d c_inv = Vector3d(1.0) - c;

	Lattice &t = *this;
	int i = cell.x, j = cell.y, k = cell.z;

	t(i, j, k) += c_inv.x * c_inv.y * c_inv.z * value;
	t(i, j, k + 1) += c_inv.x * c_inv.y * c.z * value;
	t(i, j + 1, k) += c_inv.x * c.y * c_inv.z * value;
	t(i, j + 1, k + 1) += c_inv.x * c.y * c.z * value;
	t(i + 1, j, k) += c.x * c_inv.y * c_inv.z * value;
	t(i + 1, j, k + 1) += c.x * c_inv.y * c.z * value;
	t(i + 1, j + 1, k) += c.x * c.y * c_inv.z * value;
	t(i + 1, j + 1, k + 1) += c.x * c.y * c.z * value;
	

	/*data[k1]     += c_inv.x * c_inv.y * c_inv.z * value;
	data[k2]     += c_inv.x * c_inv.y * c.z * value;
	data[k3]     += c_inv.x * c.y * c_inv.z * value;
	data[k4]     += c_inv.x * c.y * c.z * value;
	data[k1 + 1] += c.x * c_inv.y * c_inv.z * value;
	data[k2 + 1] += c.x * c_inv.y * c.z * value;
	data[k3 + 1] += c.x * c.y * c_inv.z * value;
	data[k4 + 1] += c.x * c.y * c.z * value;*/
}

YeeGrid::YeeGrid(const Vector3d &vmin, const Vector3d &vmax, const Vector3i &numInnerCells) :
	vmin(vmin),
	vmax(vmax),
	numCells(numInnerCells + Vector3i(2)),
	cellSize((vmax - vmin) / numInnerCells),
	invCellVolume(1.0 / (cellSize.x * cellSize.y * cellSize.z)),

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
	Vector3d j = pt.factor * pt.charge * pt.Velocity() * invCellVolume;

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
	double kj = 4 * M_PI * dt;
	Vector3d d = Vector3d(1.0) / cellSize;

	// Ex
	for (int i = 1; i < Ex.GetSize().x - 1; i++)
		for (int j = 1; j < Ex.GetSize().y - 1; j++)
			for (int k = 1; k < Ex.GetSize().z - 1; k++)
			{
				Ex(i, j, k) += cdt *
					((Bz(i, j, k) - Bz(i, j - 1, k)) * d.y -
					(By(i, j, k) - By(i, j, k - 1)) * d.z) - 
					kj * Jx(i, j, k);
			}

	// Ey
	for (int i = 1; i < Ey.GetSize().x - 1; i++)
		for (int j = 1; j < Ey.GetSize().y - 1; j++)
			for (int k = 1; k < Ey.GetSize().z - 1; k++)
			{
				Ey(i, j, k) -= cdt *
					((Bz(i, j, k) - Bz(i - 1, j, k)) * d.x -
					(Bx(i, j, k) - Bx(i, j, k - 1)) * d.z) -
					kj * Jy(i, j, k);
			}

	// Ez
	for (int i = 1; i < Ez.GetSize().x - 1; i++)
		for (int j = 1; j < Ez.GetSize().y - 1; j++)
			for (int k = 1; k < Ez.GetSize().z - 1; k++)
			{
				Ez(i, j, k) += cdt *
					((By(i, j, k) - By(i - 1, j, k)) * d.x -
					(Bx(i, j, k) - Bx(i, j - 1, k)) * d.y) -
					kj * Jz(i, j, k);
			}

	pbc_E();

	// Bx
	for (int i = 1; i < Bx.GetSize().x - 1; i++)
		for (int j = 1; j < Bx.GetSize().y - 1; j++)
			for (int k = 1; k < Bx.GetSize().z - 1; k++)
			{
				Bx(i, j, k) -= cdt *
					((Ez(i, j + 1, k) - Ez(i, j, k)) * d.y -
					(Ey(i, j, k + 1) - Ey(i, j, k)) * d.z);
			}

	// By
	for (int i = 1; i < By.GetSize().x - 1; i++)
		for (int j = 1; j < By.GetSize().y - 1; j++)
			for (int k = 1; k < By.GetSize().z - 1; k++)
			{
				By(i, j, k) += cdt *
					((Ez(i + 1, j, k) - Ez(i, j, k)) * d.x -
					(Ex(i, j, k + 1) - Ex(i, j, k)) * d.z);
			}

	// Bz
	for (int i = 1; i < Bz.GetSize().x - 1; i++)
		for (int j = 1; j < Bz.GetSize().y - 1; j++)
			for (int k = 1; k < Bz.GetSize().z - 1; k++)
			{
				Bz(i, j, k) -= cdt *
					((Ey(i + 1, j, k) - Ey(i, j, k)) * d.x -
					(Ex(i, j + 1, k) - Ex(i, j, k)) * d.y);
			}

	pbc_B();
}

#define BOUNDARY_LOOP(a1, a2, latt) \
	for (int a1 = 0; a1 < latt.GetSize().a1; a1++) \
	for (int a2 = 0; a2 < latt.GetSize().a2; a2++)

void YeeGrid::pbc(Lattice &l)
{
	Vector3i s = l.GetSize() - Vector3i(1);
	BOUNDARY_LOOP(x, y, l)
	{
		l(x, y, 0) = l(x, y, s.z - 1);
		l(x, y, s.z) = l(x, y, 1);
	}
	BOUNDARY_LOOP(x, z, l)
	{
		l(x, 0, z) = l(x, s.y - 1, z);
		l(x, s.y, z) = l(x, 1, z);
	}
	BOUNDARY_LOOP(y, z, l)
	{
		l(0, y, z) = l(s.x - 1, y, z);
		l(s.x, y, z) = l(1, y, z);
	}
}

void YeeGrid::pbc_J()
{
	// Jx
	Vector3i s = Jx.GetSize() - Vector3i(1);
	BOUNDARY_LOOP(x, y, Jx)
	{
		double a = Jx(x, y, 1) + Jx(x, y, s.z - 1);
		Jx(x, y, 1) = Jx(x, y, s.z - 1) = a;
	}
	BOUNDARY_LOOP(x, z, Jx)
	{
		double a = Jx(x, 1, z) + Jx(x, s.y - 1, z);
		Jx(x, 1, z) = Jx(x, s.y - 1, z) = a;
	}
	BOUNDARY_LOOP(y, z, Jx)
	{
		Jx(1, y, z) += Jx(s.x, y, z);
		Jx(s.x - 1, y, z) += Jx(0, y, z);
	}

	// Jy
	s = Jy.GetSize() - Vector3i(1);
	BOUNDARY_LOOP(x, y, Jy)
	{
		double a = Jy(x, y, 1) + Jy(x, y, s.z - 1);
		Jy(x, y, 1) = Jy(x, y, s.z - 1) = a;
	}
	BOUNDARY_LOOP(x, z, Jy)
	{
		Jy(x, 1, z) += Jy(x, s.y, z);
		Jy(x, s.y - 1, z) += Jy(x, 0, z);
	}
	BOUNDARY_LOOP(y, z, Jy)
	{
		double a = Jy(1, y, z) + Jy(s.x - 1, y, z);
		Jy(1, y, z) = Jy(s.x - 1, y, z) = a;
	}

	// Jz
	s = Jz.GetSize() - Vector3i(1);
	BOUNDARY_LOOP(x, y, Jz)
	{
		Jz(x, y, 1) += Jz(x, y, s.z);
		Jz(x, y, s.z - 1) += Jz(x, y, 0);
	}
	BOUNDARY_LOOP(x, z, Jz)
	{
		double a = Jz(x, 1, z) + Jz(x, s.y - 1, z);
		Jz(x, 1, z) = Jz(x, s.y - 1, z) = a;
	}
	BOUNDARY_LOOP(y, z, Jz)
	{
		double a = Jz(1, y, z) + Jz(s.x - 1, y, z);
		Jz(1, y, z) = Jz(s.x - 1, y, z) = a;
	}

}

#undef BOUNDARY_LOOP


void YeeGrid::pbc_E()
{
	pbc(Ex);
	pbc(Ey);
	pbc(Ez);
}

void YeeGrid::pbc_B()
{
	pbc(Bx);
	pbc(By);
	pbc(Bz);
}