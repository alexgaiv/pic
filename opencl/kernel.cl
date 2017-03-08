#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#include "grid.h"
#include "utils.h"

void CopyToGlobal(struct Grid *grid, int3 cId, int3 cId_g);

kernel void main(
	float3 vmin, float3 vmax, int3 numInnerCells,

	global float *Ex_g, global float *Ey_g, global float *Ez_g,
	global float *Bx_g, global float *By_g, global float *Bz_g,
	global float *Jx_g, global float *Jy_g, global float *Jz_g,

	local float *Ex, local float *Ey, local float *Ez,
	local float *Bx, local float *By, local float *Bz,
	local float *Jx, local float *Jy, local float *Jz)
{
	// ensemble size
	int3 ls = (int3)(get_local_size(0), get_local_size(1), get_local_size(2));

	// total number of cells
	int3 gs = (int3)(get_global_size(0), get_global_size(1), get_global_size(2));

	// ensemble id
	int3 eId = (int3)(get_group_id(0), get_group_id(1), get_group_id(2));

	// local cell id
	int3 cId = (int3)(get_local_id(0), get_local_id(1), get_local_id(2));

	// global cell id
	int3 cId_g = eId * ls + cId;

	struct Grid grid = createGrid(
		ls, gs, eId,
		vmin, vmax, numInnerCells,
		Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz,
		Ex_g, Ey_g, Ez_g, Bx_g, By_g, Bz_g, Jx_g, Jy_g, Jz_g);

	float3 cellSize = grid.cellSize;
	float x = vmin.x + (cId_g.x + 0.5) * cellSize.x;
	float E = 1111 * cos(2 * M_PI * x);

	int4 i = idx4(cId + (int3)1, grid.Ex.size);
	Ex[i.x] = E;
	Ex[i.y] = E;
	Ex[i.z] = E;
	Ex[i.w] = E;

	//*/
	if (cId.x == 0)
	{
		float x = vmin.x + (cId_g.x - 0.5) * cellSize.x;
		float E = 1111 * cos(2 * M_PI * x);
		int4 i = idx4(cId + (int3)(0, 1, 1), grid.Ex.size);
		Ex[i.x] = E;
		Ex[i.y] = E;
		Ex[i.z] = E;
		Ex[i.w] = E;
	}
	else if (cId.x == ls.x - 1)
	{
		float x = vmin.x + (cId_g.x + 1.5) * cellSize.x;
		float E = 1111 * cos(2 * M_PI * x);
		int4 i = idx4(cId + (int3)(2, 1, 1), grid.Ex.size);
		Ex[i.x] = E;
		Ex[i.y] = E;
		Ex[i.z] = E;
		Ex[i.w] = E;
	}

	if (cId.y == 0)
	{
		int4 i = idx4(cId + (int3)(1, 0, 1), grid.Ex.size);
		Ex[i.x] = E;
		Ex[i.z] = E;

	}
	else if (cId.y == ls.y - 1)
	{
		int4 i = idx4(cId + (int3)(1, 2, 1), grid.Ex.size);
		Ex[i.y] = E;
		Ex[i.w] = E;
	}

	if (cId.z == 0)
	{
		int4 i = idx4(cId + (int3)(1, 1, 0), grid.Ex.size);
		Ex[i.x] = E;
		Ex[i.y] = E;
	}
	else if (cId.z == ls.z - 1)
	{
		int4 i = idx4(cId + (int3)(1, 1, 2), grid.Ex.size);
		Ex[i.z] = E;
		Ex[i.w] = E;
	}
	//*/

	/*float x1 = vmin.x + cId_g.x * cellSize.x;
	float x2 = x1 + cellSize.x;
	float B1 = 1111 * cos(2 * M_PI * x1);
	float B2 = 1111 * cos(2 * M_PI * x2);

	int j = idx(cId, grid.Bx.size);
	Bx[j] = B1;
	Bx[j + 1] = B2;*/

	barrier(CLK_LOCAL_MEM_FENCE);

	//float3 p = (convert_float3(cId_g) + (float3)0.5) * cellSize;
	//struct FieldPoint f = grid_InterpolateField(&grid, p, cId);

	CopyToGlobal(&grid, cId, cId_g);
}

void CopyToGlobal(struct Grid *grid, int3 cId, int3 cId_g)
{
	int4 i1 = idx4(cId_g + (int3)1, grid->Ex.size_g);
	int4 i2 = idx4(cId + (int3)1, grid->Ex.size);
	grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
	grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
	grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
	grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];

	//*/
	if (cId_g.x == 0)
	{
		int4 i1 = idx4(cId_g + (int3)(0, 1, 1), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(0, 1, 1), grid->Ex.size);

		grid->Ex.data_g[i1.x] = grid->Ex.data[i2.x];
		grid->Ex.data_g[i1.y] = grid->Ex.data[i2.y];
		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
	else if (cId_g.x == grid->globalSize.x - 1)
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
	else if (cId_g.y == grid->globalSize.y - 1)
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
	else if (cId_g.z == grid->globalSize.z - 1)
	{
		int4 i1 = idx4(cId_g + (int3)(1, 1, 2), grid->Ex.size_g);
		int4 i2 = idx4(cId + (int3)(1, 1, 2), grid->Ex.size);

		grid->Ex.data_g[i1.z] = grid->Ex.data[i2.z];
		grid->Ex.data_g[i1.w] = grid->Ex.data[i2.w];
	}
	//*/

	/*i1 = idx4(cId_g + (int3)1, grid->Ey.size_g);
	i2 = idx4(cId, grid->Ey.size);
	grid->Ey.data_g[i1.x]     = grid->Ey.data[i2.x];
	grid->Ey.data_g[i1.x + 1] = grid->Ey.data[i2.x + 1];
	grid->Ey.data_g[i1.z]     = grid->Ey.data[i2.z];
	grid->Ey.data_g[i1.z + 1] = grid->Ey.data[i2.z + 1];

	i1 = idx4(cId_g + (int3)1, grid->Ez.size_g);
	i2 = idx4(cId, grid->Ez.size);
	grid->Ez.data_g[i1.x]     = grid->Ez.data[i2.x];
	grid->Ez.data_g[i1.x + 1] = grid->Ez.data[i2.x + 1];
	grid->Ez.data_g[i1.y]     = grid->Ez.data[i2.y];
	grid->Ez.data_g[i1.y + 1] = grid->Ez.data[i2.y + 1];

	int j1 = idx(cId_g + (int3)1, grid->Bx.size_g);
	int j2 = idx(cId, grid->Bx.size);
	grid->Bx.data_g[j1]     = grid->Bx.data[j2];
	grid->Bx.data_g[j1 + 1] = grid->Bx.data[j2 + 1];

	j1 = idx(cId_g + (int3)1, grid->By.size_g);
	j2 = idx(cId, grid->By.size);
	grid->By.data_g[j1]                     = grid->By.data[j2];
	grid->By.data_g[j1 + grid->By.size_g.x] = grid->By.data[j2 + grid->By.size.x];

	j1 = idx(cId_g + (int3)1, grid->Bz.size_g);
	j2 = idx(cId, grid->Bz.size);
	grid->Bz.data_g[j1] = grid->Bz.data[j2];
	grid->Bz.data_g[j1 + grid->Bz.size_g.x * grid->Bz.size_g.y] =
		grid->Bz.data[j2 + grid->Bz.size.x * grid->Bz.size.y];
		*/
}