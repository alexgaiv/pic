#ifndef _TESTS_H_
#define _TESTS_H_

#include "utils.h"
#include "grid.h"

struct ErrorStruct
{
	float abs;
	float rel;
};

struct CoordsMomentumError
{
	struct ErrorStruct r;
	struct ErrorStruct p;
};

void SetEx(struct Grid *grid, float val)
{
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
	else if (cId.x == grid->wi.local_size.x - 1)
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
	else if (cId.y == grid->wi.local_size.y - 1)
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
	else if (cId.z == grid->wi.local_size.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Ex.size);
		Ex[i.z] = val;
		Ex[i.w] = val;
	}
}

void SetBz(struct Grid *grid, float val)
{
	int3 cId = grid->wi.cell_id;
	local float *Bz = grid->Bz.data;

	int j = idx(cId + (int3)1, grid->Bz.size);
	int size_xy = grid->Bz.size.x * grid->Bz.size.y;
	Bz[j] = val;
	Bz[j + size_xy] = val;

	if (cId.x == 0)
	{
		int j = idx(cId + (int3)(0, 1, 1), grid->Bz.size);
		Bz[j] = val;
		Bz[j + size_xy] = val;
	}
	else if (cId.x == grid->wi.local_size.x - 1)
	{
		int j = idx(cId + (int3)(2, 1, 1), grid->Bz.size);
		Bz[j] = val;
		Bz[j + size_xy] = val;
	}

	if (cId.y == 0)
	{
		int j = idx(cId + (int3)(1, 0, 1), grid->Bz.size);
		Bz[j] = val;
		Bz[j + size_xy] = val;
	}
	else if (cId.y == grid->wi.local_size.y - 1)
	{
		int j = idx(cId + (int3)(1, 2, 1), grid->Bz.size);
		Bz[j] = val;
		Bz[j + size_xy] = val;
	}

	if (cId.z == 0)
	{
		int j = idx(cId + (int3)(1, 1, 0), grid->Bz.size);
		Bz[j] = val;
	}
	else if (cId.z == grid->wi.local_size.z - 1)
	{
		int j = idx(cId + (int3)(1, 1, 2), grid->Bz.size);
		Bz[j] = val;
	}
}

struct CoordsMomentumError TestBoris_1(struct Grid *grid, int steps, struct Particle *pt)
{
	struct CoordsMomentumError err;
	struct ErrorStruct *r_error = &err.r, *p_error = &err.p;

	struct Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;
	p.coords = (float3)0;
	p.momentum = (float3)0;

	const float mc = p.mass * c;
	const float E0 = 43.0;

	SetEx(grid, E0);

	float dt = mc / (p.charge * E0 * steps);

	MoveParticle2(&p, grid, dt, steps);

	float3 r_theor = (float3)0;
	r_theor.x = mc * c / (p.charge * E0) * (sqrt(2.0) - 1.0);
	float3 p_theor = (float3)(mc, 0, 0);

	r_error->abs = length(p.coords - r_theor);
	p_error->abs = length(p.momentum - p_theor);
	r_error->rel = r_error->abs / length(r_theor);
	p_error->rel = p_error->abs / length(p_theor);

	*pt = p;

	return err;
}

struct CoordsMomentumError TestBoris_2(struct Grid *grid, int steps, struct Particle *pt)
{
	struct CoordsMomentumError err;
	struct ErrorStruct *r_error = &err.r, *p_error = &err.p;

	const float mc = electronMass * c;
	const float p0 = 5.0;
	const float B0 = 57.0;

	struct Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;
	p.coords = (float3)0;
	p.momentum = (float3)(p0, 0, 0);
	
	SetBz(grid, B0);

	float dt = M_PI * mc / (fabs(p.charge) * B0 * steps) * sqrt(1 + pow(p0 / mc, 2));

	MoveParticle2(&p, grid, dt, steps);

	float3 r_theor = (float3)(0, -2 * p0 * c / (p.charge * B0), 0);
	float3 p_theor = (float3)(-p0, 0, 0);

	r_error->abs = length(p.coords - r_theor);
	p_error->abs = length(p.momentum - p_theor);
	r_error->rel = r_error->abs / length(r_theor);
	p_error->rel = p_error->abs / length(p_theor);

	*pt = p;

	return err;
}

void TestGrid(struct Grid *grid, global int *result)
{
	float3 cellSize = grid->cellSize;
	float3 vmin = grid->vmin;
	int3 cId = grid->wi.cell_id;
	int3 cId_g = grid->wi.global_cell_id;
	local float *Ex = grid->Ex.data;

	float x = vmin.x + (cId_g.x + 0.5) * cellSize.x;
	int4 i = idx4(cId + (int3)1, grid->Ex.size);
	Ex[i.x] = x;
	Ex[i.y] = x;
	Ex[i.z] = x;
	Ex[i.w] = x;

	if (cId.x == 0)
	{
		float x = vmin.x + (cId_g.x - 0.5) * cellSize.x;
		int4 i = idx4(cId + (int3)(0, 1, 1), grid->Ex.size);
		Ex[i.x] = x;
		Ex[i.y] = x;
		Ex[i.z] = x;
		Ex[i.w] = x;
	}
	else if (cId.x == grid->wi.local_size.x - 1)
	{
		float x = vmin.x + (cId_g.x + 1.5) * cellSize.x;
		int4 i = idx4(cId + (int3)(2, 1, 1), grid->Ex.size);
		Ex[i.x] = x;
		Ex[i.y] = x;
		Ex[i.z] = x;
		Ex[i.w] = x;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid->Ex.size);
		Ex[i.x] = x;
		Ex[i.z] = x;
	}
	else if (cId.y == grid->wi.local_size.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid->Ex.size);
		Ex[i.y] = x;
		Ex[i.w] = x;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid->Ex.size);
		Ex[i.x] = x;
		Ex[i.y] = x;
	}
	else if (cId.z == grid->wi.local_size.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid->Ex.size);
		Ex[i.z] = x;
		Ex[i.w] = x;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	float3 p = (convert_float3(cId_g) + (float3)0.3) * cellSize;
	struct FieldPoint f = grid_InterpolateField(grid, p);
	result[idx(cId_g, grid->wi.global_size)] = fabs(f.E.x - p.x) < PRECISION ? 1 : 0;
}

#endif // _TESTS_H_