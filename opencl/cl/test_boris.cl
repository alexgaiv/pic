#include "utils.hcl"
#include "grid.hcl"
#include "tests.hcl"

kernel void main(
    float3 vmin, float3 vmax,

    global float *Ex_g, global float *Ey_g, global float *Ez_g,
    global float *Bx_g, global float *By_g, global float *Bz_g,
    global float *Jx_g, global float *Jy_g, global float *Jz_g,

    local float *Ex, local float *Ey, local float *Ez,
    local float *Bx, local float *By, local float *Bz,
    local float *Jx, local float *Jy, local float *Jz,

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
    WorkItemInfo wi;
    initWorkItemInfo(&wi);

    Grid grid;
    initGrid(
        &grid, &wi, vmin, vmax,
        Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz,
        Ex_g, Ey_g, Ez_g, Bx_g, By_g, Bz_g, Jx_g, Jy_g, Jz_g,
        0, 0, 0, 0, 0, 0);

    TestBoris(&grid, startN, dN, numIterations,
        r_rel_test1, r_abs_test1,
        p_rel_test1, p_abs_test1,
        r_rel_test2, r_abs_test2,
        p_rel_test2, p_abs_test2);
}