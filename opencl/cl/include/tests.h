#ifndef _TESTS_H_
#define _TESTS_H_

#include "utils.h"
#include "grid.h"

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