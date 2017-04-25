#ifndef _TESTS_H_
#define _TESTS_H_

#include "utils.hcl"
#include "particle.hcl"
#include "particle_mover.hcl"
#include "grid.hcl"
#include "init_fields.hcl"

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

struct CoordsMomentumError TestBoris_1(struct Grid *grid, int steps, float E0)
{
    struct CoordsMomentumError err;
    struct ErrorStruct *r_error = &err.r, *p_error = &err.p;

    struct Particle p;
    p.mass = electronMass;
    p.charge = electronCharge;
    p.coords = (float3)0;
    particle_SetMomentum(&p, (float3)0);

    const float mc = p.mass * c;
    float dt = mc / (p.charge * E0 * steps);

    MoveParticle2(&p, grid, dt, steps);

    float3 r_theor = (float3)0;
    r_theor.x = mc * c / (p.charge * E0) * (sqrt(2.0) - 1.0);
    float3 p_theor = (float3)(mc, 0, 0);

    r_error->abs = length(p.coords - r_theor);
    p_error->abs = length(particle_GetMomentum(&p) - p_theor);
    r_error->rel = r_error->abs / length(r_theor);
    p_error->rel = p_error->abs / length(p_theor);

    return err;
}

struct CoordsMomentumError TestBoris_2(struct Grid *grid, int steps, float B0)
{
    struct CoordsMomentumError err;
    struct ErrorStruct *r_error = &err.r, *p_error = &err.p;

    const float mc = electronMass * c;
    const float p0 = 5.0;

    struct Particle p;
    p.mass = electronMass;
    p.charge = electronCharge;
    p.coords = (float3)0.0;
    particle_SetMomentum(&p, (float3)(p0, 0.0, 0.0));
    
    float dt = M_PI * mc / (fabs(p.charge) * B0 * (float)steps) * sqrt(1.0 + pow(p0 / mc, 2.0f));

    MoveParticle2(&p, grid, dt, steps);

    float3 r_theor = (float3)(0.0, -2.0 * p0 * c / (p.charge * B0), 0.0);
    float3 p_theor = (float3)(-p0, 0.0, 0.0);

    r_error->abs = length(p.coords - r_theor);
    p_error->abs = length(particle_GetMomentum(&p) - p_theor);
    r_error->rel = r_error->abs / length(r_theor);
    p_error->rel = p_error->abs / length(p_theor);

    return err;
}

void TestBoris(
    struct Grid *grid,
    int startN, int dN, int numIterations,
    global float *r_rel_test1,
    global float *r_abs_test1,
    global float *p_rel_test1,
    global float *p_abs_test1,

    global float *r_rel_test2,
    global float *r_abs_test2,
    global float *p_rel_test2,
    global float *p_abs_test2)
{
    const float E0 = 43.0;
    const float B0 = 57.0;

    SetEx(grid, E0);
    SetEy(grid, 0.0);
    SetEz(grid, 0.0);
    SetBx(grid, E0);
    SetBy(grid, 0.0);
    SetBz(grid, 0.0);
    barrier(CLK_LOCAL_MEM_FENCE);

    int3 cId_g = grid->wi.global_cell_id;
    if (cId_g.x == 0 && cId_g.y == 0 && cId_g.z == 0)
    {
        int N = startN;
        for (int i = 0; i < numIterations; i++, N += dN)
        {
            struct CoordsMomentumError err = TestBoris_1(grid, N, E0);
            r_rel_test1[i] = err.r.rel;
            r_abs_test1[i] = err.r.abs;
            p_rel_test1[i] = err.p.rel;
            p_abs_test1[i] = err.p.abs;
        }
    }

    SetEx(grid, 0);
    SetBz(grid, B0);
    barrier(CLK_LOCAL_MEM_FENCE);

    if (cId_g.x == 0 && cId_g.y == 0 && cId_g.z == 0)
    {
        int N = startN;
        for (int i = 0; i < numIterations; i++, N += dN)
        {
            struct CoordsMomentumError err = TestBoris_2(grid, N, B0);
            r_rel_test2[i] = err.r.rel;
            r_abs_test2[i] = err.r.abs;
            p_rel_test2[i] = err.p.rel;
            p_abs_test2[i] = err.p.abs;
        }
    }
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