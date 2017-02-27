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

	float3 cellSize = (vmax - vmin) / (float3)(numInnerCells.x, numInnerCells.y, numInnerCells.z);

	// arrays sizes (z is not used)
	int2 sEx = ls.xy + (int2)(0, 1);
	int2 sEy = ls.xy + (int2)(1, 0);
	int2 sEz = ls.xy + (int2)(1, 1);
	int2 sBx = ls.xy + (int2)(1, 0);
	int2 sBy = ls.xy + (int2)(0, 1);
	int2 sBz = ls.xy;

	int2 sEx_g = gs.xy + (int2)(2, 3); // + boundary cells
	int2 sEy_g = gs.xy + (int2)(3, 2);
	int2 sEz_g = gs.xy + (int2)(3, 3);
	int2 sBx_g = gs.xy + (int2)(3, 2);
	int2 sBy_g = gs.xy + (int2)(2, 3);
	int2 sBz_g = gs.xy + (int2)(2, 2);

	float x = vmin.x + (cId_g.x - 0.5) * cellSize.x;
	float E = 1111 * cos(2 * M_PI * x);

	int i1 = (cId.z * sEx.y + cId.y) * sEx.x + cId.x;
	int i2 = (cId.z * sEx.y + cId.y + 1) * sEx.x + cId.x;
	int i3 = ((cId.z + 1) * sEx.y + cId.y) * sEx.x + cId.x;
	int i4 = ((cId.z + 1) * sEx.y + cId.y + 1) * sEx.x + cId.x;

	Ex[i1] = E;
	Ex[i2] = E;
	Ex[i3] = E;
	Ex[i4] = E;

	

	barrier(CLK_LOCAL_MEM_FENCE);

	// copy local memory to global
	int3 t = cId_g + (int3)(1, 1, 1);
	Ex_g[(t.z * sEx_g.y + t.y) * sEx_g.x + t.x]           = Ex[i1];
	Ex_g[(t.z * sEx_g.y + t.y + 1) * sEx_g.x + t.x]       = Ex[i2];
	Ex_g[((t.z + 1) * sEx_g.y + t.y) * sEx_g.x + t.x]     = Ex[i3];
	Ex_g[((t.z + 1) * sEx_g.y + t.y + 1) * sEx_g.x + t.x] = Ex[i4];
}