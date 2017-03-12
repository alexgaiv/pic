#ifndef _GRID_H_
#define _GRID_H_

#include "utils.h"

constant float3 shift_EJx = (float3)(0.5, 1, 1);
constant float3 shift_EJy = (float3)(1, 0.5, 1);
constant float3 shift_EJz = (float3)(1, 1, 0.5);

constant float3 shift_Bx = (float3)(1, 0.5, 0.5);
constant float3 shift_By = (float3)(0.5, 1, 0.5);
constant float3 shift_Bz = (float3)(0.5, 0.5, 1);

struct FieldPoint
{
	float3 E;
	float3 B;
};

struct Field
{
	local float *data;
	int3 size;

	global float *data_g;
	int3 size_g;
};

void initField(
 	struct Field *field,
	local float *data, int3 size,
	global float *data_g, int3 size_g)
{
	field->data = data;
	field->size = size;
	field->data_g = data_g;
	field->size_g = size_g;
}

double field_Interpolate(struct Field *field, int3 cell, float3 coords)
{
	int4 i = idx4(cell, field->size);

	float3 c = coords - convert_float3(cell);
	float3 c_inv = (float3)1.0 - c;

	local float *data = field->data;
	double c00 = data[i.x] * c_inv.x + data[i.x + 1] * c.x;
	double c01 = data[i.y] * c_inv.x + data[i.y + 1] * c.x;
	double c10 = data[i.z] * c_inv.x + data[i.z + 1] * c.x;
	double c11 = data[i.w] * c_inv.x + data[i.w + 1] * c.x;

	double c0 = c00 * c_inv.y + c10 * c.y;
	double c1 = c01 * c_inv.y + c11 * c.y;

	return c0 * c_inv.z + c1 * c.z;
}

struct Grid
{
	struct WorkItemInfo wi;
	float3 vmin, vmax;
	float3 cellSize;

	struct Field Ex, Ey, Ez;
	struct Field Bx, By, Bz;
	struct Field Jx, Jy, Jz;
};

void initGrid(
 	struct Grid *grid,
	struct WorkItemInfo *wi,
	float3 vmin, float3 vmax,
	int3 numInnerCells,
	local float *Ex, local float *Ey, local float *Ez,
	local float *Bx, local float *By, local float *Bz,
	local float *Jx, local float *Jy, local float *Jz,

	global float *Ex_g, global float *Ey_g, global float *Ez_g,
	global float *Bx_g, global float *By_g, global float *Bz_g,
	global float *Jx_g, global float *Jy_g, global float *Jz_g)
{
	grid->vmin = vmin;
	grid->vmax = vmax;
	grid->cellSize = (vmax - vmin) / convert_float3(numInnerCells);
	grid->wi = *wi;
	
	int3 ls = wi->local_size;
	int3 gs = wi->global_size;

	int3 sEx = ls + (int3)(2, 3, 3);
	int3 sEy = ls + (int3)(3, 2, 3);
	int3 sEz = ls + (int3)(3, 3, 2);
	int3 sBx = ls + (int3)(3, 2, 2);
	int3 sBy = ls + (int3)(2, 3, 2);
	int3 sBz = ls + (int3)(2, 2, 3);

	int3 sEx_g = gs + (int3)(2, 3, 3);
	int3 sEy_g = gs + (int3)(3, 2, 3);
	int3 sEz_g = gs + (int3)(3, 3, 2);
	int3 sBx_g = gs + (int3)(3, 2, 2);
	int3 sBy_g = gs + (int3)(2, 3, 2);
	int3 sBz_g = gs + (int3)(2, 2, 3);

	initField(&grid->Ex, Ex, sEx, Ex_g, sEx_g);
	initField(&grid->Ey, Ey, sEy, Ey_g, sEy_g);
	initField(&grid->Ez, Ez, sEz, Ez_g, sEz_g);
	initField(&grid->Bx, Bx, sBx, Bx_g, sBx_g);
	initField(&grid->By, By, sBy, By_g, sBy_g);
	initField(&grid->Bz, Bz, sBz, Bz_g, sBz_g);
	initField(&grid->Jz, Jz, sEz, Jz_g, sEz_g);
	initField(&grid->Jy, Jy, sEy, Jy_g, sEy_g);
	initField(&grid->Jz, Jz, sEz, Jz_g, sEz_g);
}

struct FieldPoint grid_InterpolateField(struct Grid *grid, float3 coords, int3 cell)
{
	float3 pos = (coords - grid->vmin) / grid->cellSize -
		convert_float3(grid->wi.group_id * grid->wi.local_size);
	int3 cell2 = convert_int3(floor(pos + (float3)0.5));

	int3 cell_Ex = cell + (int3)1;
	int3 cell_Ey = cell_Ex;
	int3 cell_Ez = cell_Ex;
	cell_Ex.x = cell2.x;
	cell_Ey.y = cell2.y;
	cell_Ez.z = cell2.z;

	int3 cell_Bx = cell2;
	int3 cell_By = cell2;
	int3 cell_Bz = cell2;
	cell_Bx.x = cell.x + 1;
	cell_By.y = cell.y + 1;
	cell_Bz.z = cell.z + 1;

	struct FieldPoint f = {
		(float3)(
			field_Interpolate(&grid->Ex, cell_Ex, pos + shift_EJx),
			field_Interpolate(&grid->Ey, cell_Ey, pos + shift_EJy),
			field_Interpolate(&grid->Ez, cell_Ez, pos + shift_EJz)),
		(float3)(
			field_Interpolate(&grid->Bx, cell_Bx, pos + shift_Bx),
			field_Interpolate(&grid->By, cell_By, pos + shift_By),
			field_Interpolate(&grid->Bz, cell_Bz, pos + shift_Bz))
	};
	return f;
}

#endif // _GRID_H_