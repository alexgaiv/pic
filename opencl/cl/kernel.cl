#include "utils.hcl"
#include "grid.hcl"
#include "particle_mover.hcl"
#include "utils.hcl"

void CopyBoundsJx(Grid *grid);
void ReduceJx(Grid *grid);
void ZeroJSum(Grid *grid);
void CopyEJx(Grid *grid, Field *field);
void CopyBz(Grid *grid);

kernel void main(
    float3 vmin, float3 vmax,

    global float *Ex_g, global float *Ey_g, global float *Ez_g,
    global float *Bx_g, global float *By_g, global float *Bz_g,
    global float *Jx_g, global float *Jy_g, global float *Jz_g,

    local float *Ex, local float *Ey, local float *Ez,
    local float *Bx, local float *By, local float *Bz,
    local float *Jx, local float *Jy, local float *Jz,

    local float *Jx_sum, local float *Jy_sum, local float *Jz_sum, // per-cell local buffers
    global float *Jx_bounds, global float *Jy_bounds, global float *Jz_bounds, // inter-group global buffers

    global float *Jx_sum_g, global float *Jy_sum_g, global float *Jz_sum_g, // for debugging
    global float3 *particles_g)
{
    WorkItemInfo wi;
    initWorkItemInfo(&wi);

    Grid grid;
    initGrid(
        &grid, &wi, vmin, vmax,
        Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz,
        Ex_g, Ey_g, Ez_g, Bx_g, By_g, Bz_g, Jx_g, Jy_g, Jz_g,
        Jx_sum, Jy_sum, Jz_sum,
        Jx_bounds, Jy_bounds, Jz_bounds);

    const int particlesPerCell = 100;
    int3 cell_id_g = wi.global_cell_id;
    int offset_local = 12*idx(wi.cell_id + (int3)1, wi.local_size + (int3)2);
    int offset_global = 12*idx(wi.global_cell_id + (int3)1, wi.global_size + (int3)2);
    int pt_offset = particlesPerCell * idx(cell_id_g, wi.global_size);

    for (int i = 0; i < 12; i++)
    {
        Jx_sum_g[offset_global + i] = 0;
        Jy_sum_g[offset_global + i] = 0;
        Jz_sum_g[offset_global + i] = 0;
    }

    ZeroJSum(&grid);

    barrier(CLK_LOCAL_MEM_FENCE);

    //if (cell_id_g.x == 0 && cell_id_g.y == 0 && cell_id_g.z == 0)
    {
        for (int i = 0; i < particlesPerCell; i++)
        {
            Particle p;
            p.coords = particles_g[pt_offset + i];
            p.mass = electronMass;
            p.charge = electronCharge;
            p.factor = 1.0;
            particle_SetVelocity(&p, (float3)1);

            grid_DepositCurrents(&grid, &p);
        }
    }

    CopyBoundsJx(&grid);

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

void CopyBoundsJx(Grid *grid)
{
    int3 cell_id       = grid->wi.cell_id;
    int3 cell_id_g     = grid->wi.global_cell_id;
    int3 group_id      = grid->wi.group_id;
    int3 local_size    = grid->wi.local_size;
    int3 minor_size    = grid->Jx_sum.size;

    for (int delta = 0; delta < 2; delta++)
    {
        int bound = delta == 0 ? 0 : local_size.x - 1;
        if (cell_id.x == bound)
        {
            int sum_offset = 12*idx(cell_id + (int3)1, local_size + (int3)2);
            int3 x_bounds_size = grid->Jx_bounds.x_size;
            int3 bounds_cell = (int3)(group_id.x, cell_id_g.y, cell_id_g.z) + (int3)(delta, 1, 1);
            int bounds_offset = 12*idx(bounds_cell, x_bounds_size);

            local float *jx_sum_data = grid->Jx_sum.data + sum_offset;
            global float *jx_bounds_data = grid->Jx_bounds.data + bounds_offset;

            int4 i = idx4((int3)(delta, 0, 0), minor_size) + (int4)delta;

            jx_bounds_data[i.x] = jx_sum_data[i.x];
            jx_bounds_data[i.y] = jx_sum_data[i.y];
            jx_bounds_data[i.z] = jx_sum_data[i.z];
            jx_bounds_data[i.w] = jx_sum_data[i.w];
        }
    }

    if (cell_id.y == 0)
    {
        int3 x_size = grid->Jx_bounds.x_size;
        int offset = x_size.x * x_size.y * x_size.z;
        global float *y_bounds = grid->Jx_bounds.data + offset;

        int3 y_bounds_size = grid->Jx_bounds.y_size;
        int3 bounds_cell = (int3)(cell_id_g.x, group_id.y, cell_id_g.z) + (int3)(0, 1, 1);
        int bounds_offset = 12*idx(bounds_cell, y_bounds_size);
    }
    else if (cell_id.y == local_size.y - 1)
    {

    }

    if (cell_id.z == 0)
    {

    }
    else if (cell_id.z == local_size.z - 1)
    {

    }
}

void ReduceJx(Grid *grid)
{
    int3 cell_id       = grid->wi.cell_id + (int3)1;
    int3 group_id      = grid->wi.group_id;
    int3 local_size    = grid->wi.local_size + (int3)2;
    int3 x_bounds_size = grid->Jx_bounds.x_size;
    int3 minor_size    = grid->Jx_bounds.minor_size;

    int i = cell_id.x;
    int j = cell_id.y;
    int k = cell_id.z;

    float j1 = 0.0f;
    float j2 = 0.0f;
    float j3 = 0.0f;
    float j4 = 0.0f;

    local float *jx = grid->Jx_sum.data;
    global float *jx_bounds = grid->Jx_bounds.data;
    
    int low = 0, hi = 3;

    for (int d = 0; d < 1; d++)
    {
        int bound = d == 0 ? local_size.x - 2 : 1;
        if (cell_id.x != bound) continue;

        int gx = group_id.x + (1 - d);
        //d == 0 ? hi-- : low--;
        d == 0 ? hi-- : low++;
        
        j1 +=
            jx_bounds[12 * idx_(gx, j, k, x_bounds_size) + idx_(d, 0, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j - 1, k, x_bounds_size) + idx_(d, 1, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j, k - 1, x_bounds_size) + idx_(d, 0, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j - 1, k - 1, x_bounds_size) + idx_(d, 1, 1, minor_size)];

        j2 +=
            jx_bounds[12 * idx_(gx, j, k, x_bounds_size) + idx_(d, 1, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j + 1, k, x_bounds_size) + idx_(d, 0, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j, k - 1, x_bounds_size) + idx_(d, 1, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j + 1, k - 1, x_bounds_size) + idx_(d, 0, 1, minor_size)];

        j3 +=
            jx_bounds[12 * idx_(gx, j, k, x_bounds_size) + idx_(d, 0, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j - 1, k, x_bounds_size) + idx_(d, 1, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j, k + 1, x_bounds_size) + idx_(d, 0, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j - 1, k + 1, x_bounds_size) + idx_(d, 1, 0, minor_size)];

        j4 +=
            jx_bounds[12 * idx_(gx, j, k, x_bounds_size) + idx_(d, 1, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j + 1, k, x_bounds_size) + idx_(d, 0, 1, minor_size)] +
            jx_bounds[12 * idx_(gx, j, k + 1, x_bounds_size) + idx_(d, 1, 0, minor_size)] +
            jx_bounds[12 * idx_(gx, j + 1, k + 1, x_bounds_size) + idx_(d, 0, 0, minor_size)];
    }

    for (int d = low; d < hi; d++)
    {
        int d1 = i + d - 1;
        int d2 = 2 - d;

        j1 +=
            jx[12 * idx_(d1, j, k, local_size) + idx_(d2, 0, 0, minor_size)] +
            jx[12 * idx_(d1, j - 1, k, local_size) + idx_(d2, 1, 0, minor_size)] +
            jx[12 * idx_(d1, j, k - 1, local_size) + idx_(d2, 0, 1, minor_size)] +
            jx[12 * idx_(d1, j - 1, k - 1, local_size) + idx_(d2, 1, 1, minor_size)];

        j2 +=
            jx[12 * idx_(d1, j, k, local_size) + idx_(d2, 1, 0, minor_size)] +
            jx[12 * idx_(d1, j + 1, k, local_size) + idx_(d2, 0, 0, minor_size)] +
            jx[12 * idx_(d1, j, k - 1, local_size) + idx_(d2, 1, 1, minor_size)] +
            jx[12 * idx_(d1, j + 1, k - 1, local_size) + idx_(d2, 0, 1, minor_size)];

        j3 +=
            jx[12 * idx_(d1, j, k, local_size) + idx_(d2, 0, 1, minor_size)] +
            jx[12 * idx_(d1, j - 1, k, local_size) + idx_(d2, 1, 1, minor_size)] +
            jx[12 * idx_(d1, j, k + 1, local_size) + idx_(d2, 0, 0, minor_size)] +
            jx[12 * idx_(d1, j - 1, k + 1, local_size) + idx_(d2, 1, 0, minor_size)];

        j4 +=
            jx[12 * idx_(d1, j, k, local_size) + idx_(d2, 1, 1, minor_size)] +
            jx[12 * idx_(d1, j + 1, k, local_size) + idx_(d2, 0, 1, minor_size)] +
            jx[12 * idx_(d1, j, k + 1, local_size) + idx_(d2, 1, 0, minor_size)] +
            jx[12 * idx_(d1, j + 1, k + 1, local_size) + idx_(d2, 0, 0, minor_size)];
    }
    
    int4 jx_idx = idx4(cell_id, grid->Jx.size);
    grid->Jx.data[jx_idx.x] = j1;
    grid->Jx.data[jx_idx.y] = j2;
    grid->Jx.data[jx_idx.z] = j3;
    grid->Jx.data[jx_idx.w] = j4;
}

void _ZeroJSum(
    Grid *grid,
    int coord1, int coord2, int bound1, int bound2,
    int3 s, int i1, int i2)
{
    for (int d1 = 0; d1 <= 2; d1 += 2) {
        int bound = d1 == 0 ? 0 : bound1 - 1;
        if (coord1 != bound) continue;

        for (int d2 = 0; d2 <= 2; d2 += 2) {
            int bound = d2 == 0 ? 0 : bound2 - 1;
            if (coord2 != bound) continue;

            int d_arr[3] = { 1, 1, 1 };
            d_arr[i1] = d1;
            d_arr[i2] = d2;
            int3 d = (int3)(d_arr[0], d_arr[1], d_arr[2]);

            int offset = 12*idx(grid->wi.cell_id + d, s);
            for (int i = 0; i < 12; i++) {
                grid->Jx_sum.data[offset + i] = 0;
                grid->Jy_sum.data[offset + i] = 0;
                grid->Jz_sum.data[offset + i] = 0;
            }
        }
    }
}

void ZeroJSum(Grid *grid)
{
    int3 cell_id = grid->wi.cell_id;
    int3 local_size = grid->wi.local_size;
    int3 s = local_size + (int3)2;

    int cId[3] = { cell_id.x, cell_id.y, cell_id.z };
    int ls[3] = { local_size.x, local_size.y, local_size.z };

    int offset = 12*idx(cell_id + (int3)1, s);

    for (int i = 0; i < 12; i++)
    {
        grid->Jx_sum.data[offset + i] = 0;
        grid->Jy_sum.data[offset + i] = 0;
        grid->Jz_sum.data[offset + i] = 0;
    }

    for (int c = 0; c < 3; c++)
    {
        if (cId[c] == 0) {
            int4 d = (int4)1;
            ((int *)&d)[c] = 0;
            int offset = 12*idx(cell_id + d.xyz, s);
            for (int i = 0; i < 12; i++) {
                grid->Jx_sum.data[offset + i] = 0;
                grid->Jy_sum.data[offset + i] = 0;
                grid->Jz_sum.data[offset + i] = 0;
            }
        }
        if (cId[c] == ls[c] - 1) {
            int4 d = (int4)1;
            ((int *)&d)[c] = 2;
            int offset = 12*idx(cell_id + d.xyz, s);
            for (int i = 0; i < 12; i++) {
                grid->Jx_sum.data[offset + i] = 0;
                grid->Jy_sum.data[offset + i] = 0;
                grid->Jz_sum.data[offset + i] = 0;
            }
        }
    }

    if ((cell_id.x == 0 || cell_id.x == local_size.x - 1) &&
        (cell_id.y == 0 || cell_id.y == local_size.y - 1))
    {
        _ZeroJSum(grid, cell_id.x, cell_id.y, local_size.x, local_size.y, s, 0, 1);
    }

    if ((cell_id.y == 0 || cell_id.y == local_size.y - 1) &&
        (cell_id.z == 0 || cell_id.z == local_size.z - 1))
    {
        _ZeroJSum(grid, cell_id.y, cell_id.z, local_size.y, local_size.z, s, 1, 2);
    }

    if ((cell_id.x == 0 || cell_id.x == local_size.x - 1) &&
        (cell_id.z == 0 || cell_id.z == local_size.z - 1))
    {
        _ZeroJSum(grid, cell_id.x, cell_id.z, local_size.x, local_size.z, s, 0, 2);
    }

    if ((cell_id.x == 0 || cell_id.x == local_size.x - 1) &&
        (cell_id.y == 0 || cell_id.y == local_size.y - 1) &&
        (cell_id.z == 0 || cell_id.z == local_size.z - 1))
    {
        for (int dx = 0; dx <= 2; dx += 2) {
            int bound = dx == 0 ? 0 : local_size.x - 1;
            if (cell_id.x != bound) continue;

            for (int dy = 0; dy <= 2; dy += 2) {
                int bound = dy == 0 ? 0 : local_size.y - 1;
                if (cell_id.y != bound) continue;

                for (int dz = 0; dz <= 2; dz += 2) {
                    int bound = dz == 0 ? 0 : local_size.z - 1;
                    if (cell_id.z != bound) continue;

                    int3 d = (int3)(dx, dy, dz);
                    int offset = 12*idx(cell_id + d, s);
                    for (int i = 0; i < 12; i++) {
                        grid->Jx_sum.data[offset + i] = 0;
                        grid->Jy_sum.data[offset + i] = 0;
                        grid->Jz_sum.data[offset + i] = 0;
                    }
                }
            }
        }
    }
}

void CopyEJx(Grid *grid, Field *field)
{
    int3 cell_id = grid->wi.cell_id;
    int3 cell_id_g = grid->wi.global_cell_id;
    int3 global_size = grid->wi.global_size;

    local float *data = field->data;
    global float *data_g = field->data_g;

    int4 i1 = idx4(cell_id_g + (int3)1, field->size_g);
    int4 i2 = idx4(cell_id + (int3)1, field->size);
    data_g[i1.x] = data[i2.x];
    data_g[i1.y] = data[i2.y];
    data_g[i1.z] = data[i2.z];
    data_g[i1.w] = data[i2.w];

    if (cell_id_g.x == 0)
    {
        int4 i1 = idx4(cell_id_g + (int3)(0, 1, 1), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(0, 1, 1), field->size);
        data_g[i1.x] = data[i2.x];
        data_g[i1.y] = data[i2.y];
        data_g[i1.z] = data[i2.z];
        data_g[i1.w] = data[i2.w];
    }
    if (cell_id_g.x == global_size.x - 1)
    {
        int4 i1 = idx4(cell_id_g + (int3)(2, 1, 1), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(2, 1, 1), field->size);
        data_g[i1.x] = data[i2.x];
        data_g[i1.y] = data[i2.y];
        data_g[i1.z] = data[i2.z];
        data_g[i1.w] = data[i2.w];
    }
    
    if (cell_id_g.y == 0)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 0, 1), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 0, 1), field->size);
        data_g[i1.x] = data[i2.x];
        data_g[i1.z] = data[i2.z];
    }
    if (cell_id_g.y == global_size.y - 1)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 2, 1), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 2, 1), field->size);
        data_g[i1.y] = data[i2.y];
        data_g[i1.w] = data[i2.w];
    }

    if (cell_id_g.z == 0)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 1, 0), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 1, 0), field->size);
        data_g[i1.x] = data[i2.x];
        data_g[i1.y] = data[i2.y];
    }
    if (cell_id_g.z == global_size.z - 1)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 1, 2), field->size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 1, 2), field->size);
        data_g[i1.z] = data[i2.z];
        data_g[i1.w] = data[i2.w];
    }
}

void CopyBz(Grid *grid)
{
    int3 cell_id = grid->wi.cell_id;
    int3 cell_id_g = grid->wi.global_cell_id;
    int3 global_size = grid->wi.global_size;

    int4 i1 = idx4(cell_id_g + (int3)1, grid->Bz.size_g);
    int4 i2 = idx4(cell_id + (int3)1, grid->Bz.size);
    grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
    grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];

    if (cell_id_g.x == 0)
    {
        int4 i1 = idx4(cell_id_g + (int3)(0, 1, 1), grid->Bz.size_g);
        int4 i2 = idx4(cell_id + (int3)(0, 1, 1), grid->Bz.size);
        grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
        grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
    }
    if (cell_id_g.x == global_size.x - 1)
    {
        int4 i1 = idx4(cell_id_g + (int3)(2, 1, 1), grid->Bz.size_g);
        int4 i2 = idx4(cell_id + (int3)(2, 1, 1), grid->Bz.size);
        grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
        grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
    }

    if (cell_id_g.y == 0)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 0, 1), grid->Bz.size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 0, 1), grid->Bz.size);
        grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
        grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
    }
    if (cell_id_g.y == global_size.y - 1)
    {
        int4 i1 = idx4(cell_id_g + (int3)(1, 2, 1), grid->Bz.size_g);
        int4 i2 = idx4(cell_id + (int3)(1, 2, 1), grid->Bz.size);
        grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
        grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
    }

    if (cell_id_g.z == 0)
    {
        int i1 = idx(cell_id_g + (int3)(1, 1, 0), grid->Bz.size_g);
        int i2 = idx(cell_id + (int3)(1, 1, 0), grid->Bz.size);
        grid->Bz.data_g[i1] = grid->Bz.data[i2];
    }
    if (cell_id_g.z == global_size.z - 1)
    {
        int i1 = idx(cell_id_g + (int3)(1, 1, 2), grid->Bz.size_g);
        int i2 = idx(cell_id + (int3)(1, 1, 2), grid->Bz.size);
        grid->Bz.data_g[i1] = grid->Bz.data[i2];
    }

    if ((cell_id_g.x == 0 || cell_id_g.x == global_size.x - 1) &&
        (cell_id_g.y == 0 || cell_id_g.y == global_size.y - 1))
    {
        int dx = cell_id_g.x == 0 ? 0 : 2;
        int dy = cell_id_g.y == 0 ? 0 : 2;
        int3 delta = (int3)(dx, dy, 1);

        int4 i1 = idx4(cell_id_g + delta, grid->Bz.size_g);
        int4 i2 = idx4(cell_id + delta, grid->Bz.size);
        grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
        grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];

        /*for (int dx = 0; dx <= 2; dx += 2) {
            int bound = dx == 0 ? 0 : global_size.x - 1;
            if (cell_id_g.x != bound) continue;

            for (int dy = 0; dy <= 2; dy += 2) {
                int bound = dy == 0 ? 0 : global_size.y - 1;
                if (cell_id_g.y != bound) continue;

                int3 delta = (int3)(dx, dy, 1);
                int4 i1 = idx4(cell_id_g + delta, grid->Bz.size_g);
                int4 i2 = idx4(cell_id + delta, grid->Bz.size);
                grid->Bz.data_g[i1.x] = grid->Bz.data[i2.x];
                grid->Bz.data_g[i1.z] = grid->Bz.data[i2.z];
            }
        }*/
    }
}