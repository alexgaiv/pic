#include "utils.h"
#include "grid.h"
#include "particle_mover.h"
#include "tests.h"

void CopyEx(struct Grid *grid);
void CopyBz(struct Grid *grid);

kernel void main(
	float3 vmin, float3 vmax, int3 numInnerCells,

	global float *Ex_g, global float *Ey_g, global float *Ez_g,
	global float *Bx_g, global float *By_g, global float *Bz_g,
	global float *Jx_g, global float *Jy_g, global float *Jz_g,

	local float *Ex, local float *Ey, local float *Ez,
	local float *Bx, local float *By, local float *Bz,
	local float *Jx, local float *Jy, local float *Jz,
	global float2 *err,
	global float3 *rp)
{
	struct WorkItemInfo wi;
	initWorkItemInfo(&wi);

	struct Grid grid;
	initGrid(
		&grid, &wi, vmin, vmax, numInnerCells,
		Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz,
		Ex_g, Ey_g, Ez_g, Bx_g, By_g, Bz_g, Jx_g, Jy_g, Jz_g);


	int3 cId = wi.global_cell_id;
	if (cId.x == 0 && cId.y == 0 && cId.z == 0)
	{
		struct Particle p;
		struct CoordsMomentumError e = TestBoris_2(&grid, 100, &p);
		err[0].x = e.r.abs;
		err[0].y = e.r.rel;
		err[1].x = e.p.abs;
		err[1].y = e.p.rel;

		rp[0] = p.coords;
		rp[1] = p.momentum;
	}

	//barrier(CLK_LOCAL_MEM_FENCE);
	
	//CopyEx(&grid);
	//CopyBz(&grid);
}

void CopyEx(struct Grid *grid)
{
	int3 cId = grid->wi.cell_id;
	int3 cId_g = grid->wi.global_cell_id;

	int4 i1 = idx4(cId_g + (int3)1, grid->Ex.size_g);
	int4 i2 = idx4(cId + (int3)1, grid->Ex.size);
	grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
	grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
	grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
	grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];

	if (cId_g.x == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(0, 1, 1), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(0, 1, 1), grid->Ex.size);
		grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
		grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
	else if (cId_g.x == grid->wi.global_size.x - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(2, 1, 1), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(2, 1, 1), grid->Ex.size);
		grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
		grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
	
	if (cId_g.y == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 0, 1), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(1, 0, 1), grid->Ex.size);
		grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
	}
	else if (cId_g.y == grid->wi.global_size.y - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 2, 1), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(1, 2, 1), grid->Ex.size);
		grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
	if (cId_g.z == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 1, 0), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(1, 1, 0), grid->Ex.size);
		grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
		grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
	}
	else if (cId_g.z == grid->wi.global_size.z - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 1, 2), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(1, 1, 2), grid->Ex.size);
		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
}

void CopyBz(struct Grid *grid)
{
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
	else if (cId_g.x == grid->wi.global_size.x - 1)
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
	else if (cId_g.y == grid->wi.global_size.y - 1)
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
	else if (cId_g.z == grid->wi.global_size.z - 1)
	{
		int i1 = idx(cId_g + (int3)(1, 1, 2), grid->Bz.size_g);
		int i2 = idx(cId + (int3)(1, 1, 2), grid->Bz.size);
		grid->Bz.data_g[i1] = grid->Bz.data[i2];
	}
}