#include "utils.h"
#include "grid.h"
#include "particle_mover.h"
#include "utils.h"

void ReduceJx(struct Grid *grid);
void ZeroJ(struct Grid *grid);
void CopyEJx(struct Grid *grid, struct Field *field);
void CopyBz(struct Grid *grid);

kernel void main(
	float3 vmin, float3 vmax, int3 numInnerCells,

	global float *Ex_g, global float *Ey_g, global float *Ez_g,
	global float *Bx_g, global float *By_g, global float *Bz_g,
	global float *Jx_g, global float *Jy_g, global float *Jz_g,

	local float *Ex, local float *Ey, local float *Ez,
	local float *Bx, local float *By, local float *Bz,
	local float *Jx, local float *Jy, local float *Jz,

	local float *Jx_sum, local float *Jy_sum, local float *Jz_sum,
	global float *Jx_sum_g, global float *Jy_sum_g, global float *Jz_sum_g,

	global float3 *particles_g,

	uint2 seed_g)
{
	struct WorkItemInfo wi;
	initWorkItemInfo(&wi);

	uint2 seed = seed_g + idx(wi.global_cell_id, wi.global_size);

	struct Grid grid;
	initGrid(
		&grid, &wi, vmin, vmax, numInnerCells,
		Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz,
		Ex_g, Ey_g, Ez_g, Bx_g, By_g, Bz_g, Jx_g, Jy_g, Jz_g,
		Jx_sum, Jy_sum, Jz_sum);

	int offset_local = 12*idx(wi.cell_id + (int3)1, wi.local_size + (int3)2);
	int offset_global = 12*idx(wi.global_cell_id + (int3)1, wi.global_size + (int3)2);
	int3 cId_g = wi.global_cell_id;
	int3 cId = grid.wi.cell_id;

	for (int i = 0; i < 12; i++)
	{
		Jx_sum_g[offset_global + i] = 0;
		Jy_sum_g[offset_global + i] = 0;
		Jz_sum_g[offset_global + i] = 0;
	}

	ZeroJ(&grid);

	barrier(CLK_LOCAL_MEM_FENCE);

	int pt_offset = 100 * idx(cId_g, wi.global_size);
	
	if (cId_g.x == 0 && cId_g.y == 1 && cId_g.z == 0)
	{
		for (int i = 0; i < 1; i++)
		{
			struct Particle p;
			p.coords = particles_g[pt_offset + i];
			p.mass = electronMass;
			p.charge = electronCharge;
			p.factor = 1.0;
			particle_SetVelocity(&p, (float3)1);

			grid_DepositCurrents(&grid, &p);
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	for (int i = 0; i < 12; i++)
	{
		Jx_sum_g[offset_global + i] = Jx_sum[offset_local + i];
		Jy_sum_g[offset_global + i] = Jy_sum[offset_local + i];
		Jz_sum_g[offset_global + i] = Jz_sum[offset_local + i];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	ReduceJx(&grid);

	barrier(CLK_LOCAL_MEM_FENCE);

	CopyEJx(&grid, &grid.Jx);
}

void ReduceJx(struct Grid *grid)
{
	int3 cId = grid->wi.cell_id + (int3)1;
	int3 ls = grid->wi.local_size + (int3)2;

	int i = cId.x;
	int j = cId.y;
	int k = cId.z;

	float j1 = 0.0f;
	float j2 = 0.0f;
	float j3 = 0.0f;
	float j4 = 0.0f;

	local float *jx = grid->Jx_sum.data;
	int3 jxSize = (int3)(3, 2, 2);

	for (int d = 0; d < 3; d++)
	{
		int d1 = i + d - 1;
		int d2 = 2 - d;

		j1 +=
			jx[12 * idx_(d1, j, k, ls) + idx_(d2, 0, 0, jxSize)] +
			jx[12 * idx_(d1, j - 1, k, ls) + idx_(d2, 1, 0, jxSize)] +
			jx[12 * idx_(d1, j, k - 1, ls) + idx_(d2, 0, 1, jxSize)] +
			jx[12 * idx_(d1, j - 1, k - 1, ls) + idx_(d2, 1, 1, jxSize)];

		j2 +=
			jx[12 * idx_(d1, j, k, ls) + idx_(d2, 1, 0, jxSize)] +
			jx[12 * idx_(d1, j + 1, k, ls) + idx_(d2, 0, 0, jxSize)] +
			jx[12 * idx_(d1, j, k - 1, ls) + idx_(d2, 1, 1, jxSize)] +
			jx[12 * idx_(d1, j + 1, k - 1, ls) + idx_(d2, 0, 1, jxSize)];

		j3 +=
			jx[12 * idx_(d1, j, k, ls) + idx_(d2, 0, 1, jxSize)] +
			jx[12 * idx_(d1, j - 1, k, ls) + idx_(d2, 1, 1, jxSize)] +
			jx[12 * idx_(d1, j, k + 1, ls) + idx_(d2, 0, 0, jxSize)] +
			jx[12 * idx_(d1, j - 1, k + 1, ls) + idx_(d2, 1, 0, jxSize)];

		j4 +=
			jx[12 * idx_(d1, j, k, ls) + idx_(d2, 1, 1, jxSize)] +
			jx[12 * idx_(d1, j + 1, k, ls) + idx_(d2, 0, 1, jxSize)] +
			jx[12 * idx_(d1, j, k + 1, ls) + idx_(d2, 1, 0, jxSize)] +
			jx[12 * idx_(d1, j + 1, k + 1, ls) + idx_(d2, 0, 0, jxSize)];
	}

	int4 i1 = idx4(cId, grid->Jx.size);
	grid->Jx.data[i1.x] = j1;
	grid->Jx.data[i1.y] = j2;
	grid->Jx.data[i1.z] = j3;
	grid->Jx.data[i1.w] = j4;
}

void ZeroJ(struct Grid *grid)
{
	int3 cId = grid->wi.cell_id;
	int3 ls = grid->wi.local_size;
	int3 s = ls + (int3)2;

	int cell_id[3] = { cId.x, cId.y, cId.z };
	int local_size[3] = { ls.x, ls.y, ls.z };

	int offset = 12*idx(cId + (int3)1, s);

	for (int i = 0; i < 12; i++)
	{
		grid->Jx_sum.data[offset + i] = 0;
		grid->Jy_sum.data[offset + i] = 0;
		grid->Jz_sum.data[offset + i] = 0;
	}

	for (int c = 0; c < 3; c++)
	{
		if (cell_id[c] == 0) {
			int4 d = (int4)1;
			((int *)&d)[c] = 0;
			int offset = 12*idx(cId + d.xyz, s);
			for (int i = 0; i < 12; i++) {
				grid->Jx_sum.data[offset + i] = 0;
				grid->Jy_sum.data[offset + i] = 0;
				grid->Jz_sum.data[offset + i] = 0;
			}
		}
		else if (cell_id[c] == local_size[c] - 1) {
			int4 d = (int4)1;
			((int *)&d)[c] = 2;
			int offset = 12*idx(cId + d.xyz, s);
			for (int i = 0; i < 12; i++) {
				grid->Jx_sum.data[offset + i] = 0;
				grid->Jy_sum.data[offset + i] = 0;
				grid->Jz_sum.data[offset + i] = 0;
			}
		}
	}

	if ((cId.x == 0 || cId.x == ls.x - 1) && (cId.y == 0 || cId.y == ls.y - 1))
	{
		int dx = cId.x == 0 ? 0 : 2;
		int dy = cId.y == 0 ? 0 : 2;
		int3 d = (int3)(dx, dy, 1);

		int offset = 12*idx(cId + d, s);
		for (int i = 0; i < 12; i++) {
			grid->Jx_sum.data[offset + i] = 0;
			grid->Jy_sum.data[offset + i] = 0;
			grid->Jz_sum.data[offset + i] = 0;
		}
	}

	if ((cId.y == 0 || cId.y == ls.y - 1) && (cId.z == 0 || cId.z == ls.z - 1))
	{
		int dy = cId.y == 0 ? 0 : 2;
		int dz = cId.z == 0 ? 0 : 2;
		int3 d = (int3)(1, dy, dz);

		int offset = 12*idx(cId + d, s);
		for (int i = 0; i < 12; i++) {
			grid->Jx_sum.data[offset + i] = 0;
			grid->Jy_sum.data[offset + i] = 0;
			grid->Jz_sum.data[offset + i] = 0;
		}
	}

	if ((cId.x == 0 || cId.x == ls.x - 1) && (cId.z == 0 || cId.z == ls.z - 1))
	{
		int dx = cId.x == 0 ? 0 : 2;
		int dz = cId.z == 0 ? 0 : 2;
		int3 d = (int3)(dx, 1, dz);

		int offset = 12*idx(cId + d, s);
		for (int i = 0; i < 12; i++) {
			grid->Jx_sum.data[offset + i] = 0;
			grid->Jy_sum.data[offset + i] = 0;
			grid->Jz_sum.data[offset + i] = 0;
		}
	}

	if ((cId.x == 0 || cId.x == ls.x - 1) &&
		(cId.y == 0 || cId.y == ls.y - 1) &&
		(cId.z == 0 || cId.z == ls.z - 1))
	{
		int dx = cId.x == 0 ? 0 : 2;
		int dy = cId.y == 0 ? 0 : 2;
		int dz = cId.z == 0 ? 0 : 2;
		int3 d = (int3)(dx, dy, dz);

		int offset = 12*idx(cId + d, s);
		for (int i = 0; i < 12; i++) {
			grid->Jx_sum.data[offset + i] = 0;
			grid->Jy_sum.data[offset + i] = 0;
			grid->Jz_sum.data[offset + i] = 0;
		}
	}
}

void CopyEJx(struct Grid *grid, struct Field *field)
{
	int3 gs = grid->wi.global_size;
	int3 cId = grid->wi.cell_id;
	int3 cId_g = grid->wi.global_cell_id;

	local float *data = field->data;
	global float *data_g = field->data_g;

	int4 i1 = idx4(cId_g + (int3)1, field->size_g);
	int4 i2 = idx4(cId + (int3)1, field->size);
	data_g[i1.x] = data[i2.x];
	data_g[i1.y] = data[i2.y];
	data_g[i1.z] = data[i2.z];
	data_g[i1.w] = data[i2.w];

	if (cId_g.x == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(0, 1, 1), field->size_g);
		int4 i2 = idx4(cId + (int3)(0, 1, 1), field->size);
		data_g[i1.x] = data[i2.x];
		data_g[i1.y] = data[i2.y];
		data_g[i1.z] = data[i2.z];
		data_g[i1.w] = data[i2.w];
	}
	else if (cId_g.x == gs.x - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(2, 1, 1), field->size_g);
		int4 i2 = idx4(cId + (int3)(2, 1, 1), field->size);
		data_g[i1.x] = data[i2.x];
		data_g[i1.y] = data[i2.y];
		data_g[i1.z] = data[i2.z];
		data_g[i1.w] = data[i2.w];
	}
	
	if (cId_g.y == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 0, 1), field->size_g);
		int4 i2 = idx4(cId + (int3)(1, 0, 1), field->size);
		data_g[i1.x] = data[i2.x];
		data_g[i1.z] = data[i2.z];
	}
	else if (cId_g.y == gs.y - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 2, 1), field->size_g);
		int4 i2 = idx4(cId + (int3)(1, 2, 1), field->size);
		data_g[i1.y] = data[i2.y];
		data_g[i1.w] = data[i2.w];
	}

	if (cId_g.z == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 1, 0), field->size_g);
		int4 i2 = idx4(cId + (int3)(1, 1, 0), field->size);
		data_g[i1.x] = data[i2.x];
		data_g[i1.y] = data[i2.y];
	}
	else if (cId_g.z == gs.z - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 1, 2), field->size_g);
		int4 i2 = idx4(cId + (int3)(1, 1, 2), field->size);
		data_g[i1.z] = data[i2.z];
		data_g[i1.w] = data[i2.w];
	}
}

void CopyBz(struct Grid *grid)
{
	int3 gs = grid->wi.global_size;
	int3 cId = grid->wi.cell_id;
	int3 cId_g = grid->wi.global_cell_id;

	int4 i1 = idx4(cId_g + (int3)1, grid->Bz.size_g);
	int4 i2 = idx4(cId + (int3)1, grid->Bz.size);
	grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
	grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];

	if (cId_g.x == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(0, 1, 1), grid->Bz.size_g);
		int4 i2 = idx4(cId + (int3)(0, 1, 1), grid->Bz.size);
		grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
		grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
	}
	else if (cId_g.x == gs.x - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(2, 1, 1), grid->Bz.size_g);
		int4 i2 = idx4(cId + (int3)(2, 1, 1), grid->Bz.size);
		grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
		grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
	}

	if (cId_g.y == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 0, 1), grid->Bz.size_g);
		int4 i2 = idx4(cId + (int3)(1, 0, 1), grid->Bz.size);
		grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
		grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
	}
	else if (cId_g.y == gs.y - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 2, 1), grid->Bz.size_g);
		int4 i2 = idx4(cId + (int3)(1, 2, 1), grid->Bz.size);
		grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
		grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
	}

	if (cId_g.z == 0)
	{
		int i1 = idx(cId_g + (int3)(1, 1, 0), grid->Bz.size_g);
		int i2 = idx(cId + (int3)(1, 1, 0), grid->Bz.size);
		grid->Bz.data_g[i1] = grid->Bz.data[i2];
	}
	else if (cId_g.z == gs.z - 1)
	{
		int i1 = idx(cId_g + (int3)(1, 1, 2), grid->Bz.size_g);
		int i2 = idx(cId + (int3)(1, 1, 2), grid->Bz.size);
		grid->Bz.data_g[i1] = grid->Bz.data[i2];
	}

	if ((cId_g.x == 0 || cId_g.x == gs.x - 1) && (cId_g.y == 0 || cId_g.y == gs.y - 1))
	{
		int dx = cId_g.x == 0 ? 0 : 2;
		int dy = cId_g.y == 0 ? 0 : 2;
		int3 delta = (int3)(dx, dy, 1);

		int4 i1 = idx4(cId_g + delta, grid->Bz.size_g);
		int4 i2 = idx4(cId + delta, grid->Bz.size);
		grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
		grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
	}
}