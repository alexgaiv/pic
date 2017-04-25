#ifndef _GRID_H_
#define _GRID_H_

#include "utils.hcl"
#include "particle.hcl"

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

float field_Interpolate(struct Field *field, int3 cell, float3 coords)
{
    int4 i = idx4(cell, field->size);
    float3 c = coords - convert_float3(cell);
    float3 c_inv = (float3)1.0 - c;

    local float *data = field->data;
    float c00 = data[i.x] * c_inv.x + data[i.x + 1] * c.x;
    float c01 = data[i.y] * c_inv.x + data[i.y + 1] * c.x;
    float c10 = data[i.z] * c_inv.x + data[i.z + 1] * c.x;
    float c11 = data[i.w] * c_inv.x + data[i.w + 1] * c.x;

    float c0 = c00 * c_inv.y + c10 * c.y;
    float c1 = c01 * c_inv.y + c11 * c.y;

    return c0 * c_inv.z + c1 * c.z;
}

void field_Deposit(struct Field *field, int offset, int3 cell, float3 coords, float value)
{
    int4 i = idx4(cell, field->size);
    float3 c = coords - convert_float3(cell);
    float3 c_inv = (float3)1.0 - c;

    local float *data = field->data + offset;

    data[i.x]     += 1;//c_inv.x * c_inv.y * c_inv.z * value;
    data[i.y]     += 2;//c_inv.x * c_inv.y * c.z * value;
    data[i.z]     += 3;//c_inv.x * c.y * c_inv.z * value;
    data[i.w]     += 4;//c_inv.x * c.y * c.z * value;
    data[i.x + 1] += 5;//c.x * c_inv.y * c_inv.z * value;
    data[i.y + 1] += 6;//c.x * c_inv.y * c.z * value;
    data[i.z + 1] += 7;//c.x * c.y * c_inv.z * value;
    data[i.w + 1] += 8;//c.x * c.y * c.z * value;
}

struct Grid
{
    struct WorkItemInfo wi;
    float3 vmin, vmax;
    float3 cellSize;

    struct Field Ex, Ey, Ez;
    struct Field Bx, By, Bz;
    struct Field Jx, Jy, Jz;
    struct Field Jx_sum, Jy_sum, Jz_sum;

    float invCellVolume;
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
    global float *Jx_g, global float *Jy_g, global float *Jz_g,

    local float *Jx_sum, local float *Jy_sum, local float *Jz_sum)
{
    grid->vmin = vmin;
    grid->vmax = vmax;
    grid->cellSize = (vmax - vmin) / convert_float3(numInnerCells);
    grid->invCellVolume = 1.0 / (grid->cellSize.x * grid->cellSize.y * grid->cellSize.z);
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
    initField(&grid->Jx, Jx, sEx, Jx_g, sEx_g);
    initField(&grid->Jy, Jy, sEy, Jy_g, sEy_g);
    initField(&grid->Jz, Jz, sEz, Jz_g, sEz_g);

    initField(&grid->Jx_sum, Jx_sum, (int3)(3, 2, 2), 0, (int3)0);
    initField(&grid->Jy_sum, Jy_sum, (int3)(2, 3, 2), 0, (int3)0);
    initField(&grid->Jz_sum, Jz_sum, (int3)(2, 2, 3), 0, (int3)0);
}

struct FieldPoint grid_InterpolateField(struct Grid *grid, float3 coords)
{
    float3 pos = (coords - grid->vmin) / grid->cellSize -
        convert_float3(grid->wi.group_offset);
    int3 cell = grid->wi.cell_id;
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

void grid_DepositCurrents(struct Grid *grid, struct Particle *pt)
{
    int offset = 12*idx(grid->wi.cell_id + (int3)1, grid->wi.local_size + (int3)2);

    float3 j = pt->factor * pt->charge * particle_GetVelocity(pt) * grid->invCellVolume;
    float3 pos = (pt->coords - grid->vmin) / grid->cellSize - 
        convert_float3(grid->wi.global_cell_id);
    int3 cell2 = convert_int3(floor(pos + (float3)0.5));

    int3 cell_Jx = (int3)(cell2.x, 0, 0);
    int3 cell_Jy = (int3)(0, cell2.y, 0);
    int3 cell_Jz = (int3)(0, 0, cell2.z);

    field_Deposit(&grid->Jx_sum, offset, cell_Jx, pos + (float3)(0.5, 0.0, 0.0), j.x);
    field_Deposit(&grid->Jy_sum, offset, cell_Jy, pos + (float3)(0.0, 0.5, 0.0), j.y);
    field_Deposit(&grid->Jz_sum, offset, cell_Jz, pos + (float3)(0.0, 0.0, 0.5), j.z);
}
#endif // _GRID_H_