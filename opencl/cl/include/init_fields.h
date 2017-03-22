#ifndef _INIT_FIELDS_H_
#define _INIT_FIELDS_H_

#include "utils.h"
#include "grid.h"

void SetEx(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *Ex = grid->Ex.data;

	int4 i = idx4(cId + (int3)1, grid->Ex.size);
	Ex[i.x] = val;
	Ex[i.y] = val;
	Ex[i.z] = val;
	Ex[i.w] = val;

	if (cId.x == 0)
	{
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->Ex.size);
		Ex[i.x] = val;
		Ex[i.y] = val;
		Ex[i.z] = val;
		Ex[i.w] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->Ex.size);
		Ex[i.x] = val;
		Ex[i.y] = val;
		Ex[i.z] = val;
		Ex[i.w] = val;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Ex.size);
		Ex[i.x] = val;
		Ex[i.z] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Ex.size);
		Ex[i.y] = val;
		Ex[i.w] = val;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->Ex.size);
		Ex[i.x] = val;
		Ex[i.y] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Ex.size);
		Ex[i.z] = val;
		Ex[i.w] = val;
	}
}

void SetEy(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *Ey = grid->Ey.data;

	int4 i = idx4(cId + (int3)1, grid->Ey.size);
	Ey[i.x] = val;
	Ey[i.z] = val;
	Ey[i.x + 1] = val;
	Ey[i.z + 1] = val;

	if (cId.x == 0)
	{
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->Ey.size);
		Ey[i.x] = val;
		Ey[i.z] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->Ey.size);
		Ey[i.x + 1] = val;
		Ey[i.z + 1] = val;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Ey.size);
		Ey[i.x] = val;
		Ey[i.z] = val;
		Ey[i.x + 1] = val;
		Ey[i.z + 1] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Ey.size);
		Ey[i.x] = val;
		Ey[i.z] = val;
		Ey[i.x + 1] = val;
		Ey[i.z + 1] = val;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->Ey.size);
		Ey[i.x] = val;
		Ey[i.x + 1] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Ey.size);
		Ey[i.z] = val;
		Ey[i.z + 1] = val;
	}
}

void SetEz(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *Ez = grid->Ez.data;

	int4 i = idx4(cId + (int3)1, grid->Ez.size);
	Ez[i.x] = val;
	Ez[i.y] = val;
	Ez[i.x + 1] = val;
	Ez[i.y + 1] = val;

	if (cId.x == 0)
	{
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->Ez.size);
		Ez[i.x] = val;
		Ez[i.y] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->Ez.size);
		Ez[i.x + 1] = val;
		Ez[i.y + 1] = val;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Ez.size);
		Ez[i.x] = val;
		Ez[i.x + 1] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Ez.size);
		Ez[i.y] = val;
		Ez[i.y + 1] = val;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->Ez.size);
		Ez[i.x] = val;
		Ez[i.y] = val;
		Ez[i.x + 1] = val;
		Ez[i.y + 1] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Ez.size);
		Ez[i.x] = val;
		Ez[i.y] = val;
		Ez[i.x + 1] = val;
		Ez[i.y + 1] = val;
	}
}

void SetBx(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *Bx = grid->Bx.data;
	
	int4 i = idx4(cId + (int3)1, grid->Bx.size);
	Bx[i.x] = val;
	Bx[i.x + 1] = val;

	if (cId.x == 0)
	{
		int i = idx(cId + (int3)(0, 1, 1), grid->Bx.size);
		Bx[i] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int i = idx(cId + (int3)(2, 1, 1), grid->Bx.size);
		Bx[i] = val;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Bx.size);
		Bx[i.x] = val;
		Bx[i.x + 1] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Bx.size);
		Bx[i.x] = val;
		Bx[i.x + 1] = val;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->Bx.size);
		Bx[i.x] = val;
		Bx[i.x + 1] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Bx.size);
		Bx[i.x] = val;
		Bx[i.x + 1] = val;
	}

	if ((cId.y == 0 || cId.y == ls.y - 1) && (cId.z == 0 || cId.z == ls.z - 1))
	{
		int dy = cId.y == 0 ? 0 : 2;
		int dz = cId.z == 0 ? 0 : 2;
		int3 delta = (int3)(1, dy, dz);

		int4 i = idx4(cId + delta, grid->Bx.size);
		Bx[i.x] = val;
		Bx[i.x + 1] = val;
	}
}


void SetBy(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *By = grid->By.data;
	
	int4 i = idx4(cId + (int3)1, grid->By.size);
	By[i.x] = val;
	By[i.y] = val;

	if (cId.x == 0)
	{
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->By.size);
		By[i.x] = val;
		By[i.y] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->By.size);
		By[i.x] = val;
		By[i.y] = val;
	}

	if (cId.y == 0)
	{
		int i = idx(cId + (int3)(1, 0, 1), grid->By.size);
		By[i] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int i = idx(cId + (int3)(1, 2, 1), grid->By.size);
		By[i] = val;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->By.size);
		By[i.x] = val;
		By[i.y] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->By.size);
		By[i.x] = val;
		By[i.y] = val;
	}

	if ((cId.x == 0 || cId.x == ls.x - 1) && (cId.z == 0 || cId.z == ls.z - 1))
	{
		int dx = cId.x == 0 ? 0 : 2;
		int dz = cId.z == 0 ? 0 : 2;
		int3 delta = (int3)(dx, 1, dz);

		int4 i = idx4(cId + delta, grid->By.size);
		By[i.x] = val;
		By[i.y] = val;
	}
}

void SetBz(struct Grid *grid, float val)
{
	int3 ls = grid->wi.local_size;
	int3 cId = grid->wi.cell_id;
	local float *Bz = grid->Bz.data;
	
	int4 i = idx4(cId + (int3)1, grid->Bz.size);
	Bz[i.x] = val;
	Bz[i.z] = val;

	if (cId.x == 0)
	{
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->Bz.size);
		Bz[i.x] = val;
		Bz[i.z] = val;
	}
	else if (cId.x == ls.x - 1)
	{
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->Bz.size);
		Bz[i.x] = val;
		Bz[i.z] = val;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Bz.size);
		Bz[i.x] = val;
		Bz[i.z] = val;
	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Bz.size);
		Bz[i.x] = val;
		Bz[i.z] = val;
	}

	if (cId.z == 0)
	{
		int i = idx(cId + (int3)(1, 1, 0), grid->Bz.size);
		Bz[i] = val;
	}
	else if (cId.z == ls.z - 1)
	{
		int i = idx(cId + (int3)(1, 1, 2), grid->Bz.size);
		Bz[i] = val;
	}

	if ((cId.x == 0 || cId.x == ls.x - 1) && (cId.y == 0 || cId.y == ls.y - 1))
	{
		int dx = cId.x == 0 ? 0 : 2;
		int dy = cId.y == 0 ? 0 : 2;
		int3 delta = (int3)(dx, dy, 1);

		int4 i = idx4(cId + delta, grid->Bz.size);
		Bz[i.x] = val;
		Bz[i.z] = val;
	}
}

#endif // _INIT_FIELDS_H_