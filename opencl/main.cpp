#include "mgl_suppress_warnings.h"
#include "common.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <mgl2\mgl.h>
#include <grid.h>
#include "pic_kernel.h"
#include "tests.h"

#pragma comment(lib, "Winmm.lib")

using std::cout;
using std::endl;

inline int idx(int i, int j, int k, const Vector3i &s) {
    return (k * s.y + j) * s.x + i;
}

inline float frand() {
    float f;
    do {
        f = (float)rand() / RAND_MAX;
    } while (f == 0.0 || f == 1.0);
    return f;
}

void ReduceJx(
    const std::vector<float> &jx, const std::vector<float> &jy, const std::vector<float> &jz,
    YeeGrid &grid)
{
    Vector3i jxSize(3, 2, 2);
    Vector3i jySize(2, 3, 2);
    Vector3i jzSize(2, 2, 3);

    Vector3i nc = grid.GetNumInnerCells() + Vector3i(2);

    FOR3_(i, j, k, nc)
    {
        float j1 = 0.0f;
        float j2 = 0.0f;
        float j3 = 0.0f;
        float j4 = 0.0f;

        for (int d = 0; d < 3; d++)
        {
            int d1 = i + d - 1;
            int d2 = 2 - d;

            j1 +=
                jx[12 * idx(d1, j, k, nc) + idx(d2, 0, 0, jxSize)] +
                jx[12 * idx(d1, j - 1, k, nc) + idx(d2, 1, 0, jxSize)] +
                jx[12 * idx(d1, j, k - 1, nc) + idx(d2, 0, 1, jxSize)] +
                jx[12 * idx(d1, j - 1, k - 1, nc) + idx(d2, 1, 1, jxSize)];

            j2 +=
                jx[12 * idx(d1, j, k, nc) + idx(d2, 1, 0, jxSize)] +
                jx[12 * idx(d1, j + 1, k, nc) + idx(d2, 0, 0, jxSize)] +
                jx[12 * idx(d1, j, k - 1, nc) + idx(d2, 1, 1, jxSize)] +
                jx[12 * idx(d1, j + 1, k - 1, nc) + idx(d2, 0, 1, jxSize)];

            j3 +=
                jx[12 * idx(d1, j, k, nc) + idx(d2, 0, 1, jxSize)] +
                jx[12 * idx(d1, j - 1, k, nc) + idx(d2, 1, 1, jxSize)] +
                jx[12 * idx(d1, j, k + 1, nc) + idx(d2, 0, 0, jxSize)] +
                jx[12 * idx(d1, j - 1, k + 1, nc) + idx(d2, 1, 0, jxSize)];

            j4 +=
                jx[12 * idx(d1, j, k, nc) + idx(d2, 1, 1, jxSize)] +
                jx[12 * idx(d1, j + 1, k, nc) + idx(d2, 0, 1, jxSize)] +
                jx[12 * idx(d1, j, k + 1, nc) + idx(d2, 1, 0, jxSize)] +
                jx[12 * idx(d1, j + 1, k + 1, nc) + idx(d2, 0, 0, jxSize)];
        }

        grid.Jx(i, j, k) = j1;
        grid.Jx(i, j + 1, k) = j2;
        grid.Jx(i, j, k + 1) = j3;
        grid.Jx(i, j + 1, k + 1) = j4;
    }
}

int main()
{
    timeBeginPeriod(1);
    //srand((unsigned)time(NULL));

    try
    {
        cout << "init context...";
        cl::Context ctx(CL_DEVICE_TYPE_CPU);
        cout << "done\n";

        cl_Descriptor cld(ctx);

        //TestGrid(cld);
        //TestBoris(cld);
        //cout << "end";
        //getchar();
        //return 0;

        YeeGrid grid0(Vector3f(0), Vector3f(1), Vector3i(4, 4, 4));
        YeeGrid grid1 = grid0;
        cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(4, 4, 4), Vector3i(2, 1, 1));

        PicKernel kernel(cld, "cl/kernel.cl", grid, true);

        Vector3f cellSize = grid.GetCellSize();
        Vector3i groupSize = grid.GetGroupSize();
        Vector3i numCells = grid.GetNumInnerCells();
        Vector3i numGroups = grid.GetNumGroups();
        int totalNumCells = numCells.x * numCells.y * numCells.z;
        int groupNumCells = groupSize.x * groupSize.y * groupSize.z;

        const int particlesPerCell = 100;
        const int numParticles = particlesPerCell * totalNumCells;
        int j_localSize = 12 * (groupSize.x + 2) * (groupSize.y + 2) * (groupSize.z + 2);
        int j_globalSize = 12 * (numCells.x + 2) * (numCells.y + 2) * (numCells.z + 2);

        cl_Buffer<cl_float3> particlesBuffer(cld, CL_MEM_WRITE_ONLY, numParticles);
        cl_Buffer<float> jxBuffer(cld, CL_MEM_WRITE_ONLY, j_globalSize);
        cl_Buffer<float> jyBuffer(cld, CL_MEM_WRITE_ONLY, j_globalSize);
        cl_Buffer<float> jzBuffer(cld, CL_MEM_WRITE_ONLY, j_globalSize);

        std::vector<cl_float3> &particles = particlesBuffer.data;
        std::vector<float> &jx = jxBuffer.data;
        std::vector<float> &jy = jyBuffer.data;
        std::vector<float> &jz = jzBuffer.data;

        FOR3(i, j, k, numCells)
        {
            Vector3f ijk((float)i, (float)j, (float)k);
            int offset = particlesPerCell * idx(i, j, k, numCells);

            Vector3f v[8] = {
                Vector3f(0.1f, 0.1f, 0.1f),
                Vector3f(0.1f, 0.9f, 0.1f),
                Vector3f(0.1f, 0.1f, 0.9f),
                Vector3f(0.1f, 0.9f, 0.9f),
                Vector3f(0.9f, 0.1f, 0.1f),
                Vector3f(0.9f, 0.9f, 0.1f),
                Vector3f(0.9f, 0.1f, 0.9f),
                Vector3f(0.9f, 0.9f, 0.9f)
            };

            for (int k = 0; k < 8; k++)
                particles[offset + k] = v2v(grid.GetMin() + (ijk + v[k]) * cellSize);

            for (int p = 8; p < particlesPerCell; p++)
            {
                particles[offset + p] = v2v(grid.GetMin() +
                    (ijk + Vector3f(frand(), frand(), frand())) * cellSize);
            }
        }

        particlesBuffer.Write();

        int boundsBufferSize = 12 *
            ((numGroups.x + 1) * (numCells.y + 2) * (numCells.z + 2) +
             (numGroups.y + 1) * (numCells.x + 2) * (numCells.z + 2) +
             (numGroups.z + 1) * (numCells.x + 2) * (numCells.y + 2));

        kernel.AddLocalBufferArg<float>(j_localSize);
        kernel.AddLocalBufferArg<float>(j_localSize);
        kernel.AddLocalBufferArg<float>(j_localSize);
        kernel.AddGlobalBufferArg<float>(boundsBufferSize, CL_MEM_READ_WRITE);
        kernel.AddGlobalBufferArg<float>(1, CL_MEM_READ_WRITE);
        kernel.AddGlobalBufferArg<float>(1, CL_MEM_READ_WRITE);
        kernel.AddArg(jxBuffer);
        kernel.AddArg(jyBuffer);
        kernel.AddArg(jzBuffer);
        kernel.AddArg(particlesBuffer);

        kernel.Run();

        jxBuffer.Read();
        jyBuffer.Read();
        jzBuffer.Read();
        grid.Jx.ReadBuffer();

        Vector3i nc = numCells + Vector3i(2);
        FOR3(i, j, k, nc)
        {
            if (i == 0 || j == 0 || k == 0 ||
                i == nc.x - 1 || j == nc.y - 1 || k == nc.z - 1)
            {
                for (int p = 0; p < 12; p++)
                    jx[12 * idx(i, j, k, nc) + p] = 0.0f;
            }
        }

        ReduceJx(jx, jy, jz, grid0);

        FOR3(i, j, k, numCells)
        {
            int offset = particlesPerCell * idx(i, j, k, numCells);
            for (int c = 0; c < particlesPerCell; c++)
            {
                cl_float4 coords = particles[offset + c];
                Particle p;
                p.coords = Vector3d(coords.x, coords.y, coords.z);
                p.charge = electronCharge;
                p.mass = electronMass;
                p.factor = 1.0;

                Vector3d velocity(1.0);
                p.momentum = velocity / sqrt(c * c - velocity.Square()) * p.mass * c;

                grid1.DepositCurrents(p);
            }
        }

        FOR3_(i, j, k, grid1.Jx.GetSize())
        {
            float j1 = (float)grid0.Jx(i, j, k);
            float j2 = (float)grid1.Jx(i, j, k);

            if (j1 == 0.0 && j2 == 0.0) continue;

            double err = abs(j1 - j2) / abs(j2);
            if (err > 0.0001) __debugbreak();
        }

        cout << "end";
    }
    catch (const cl::Error &e)
    {
        cout << endl << e.what() << " failed, err=" << e.err();
    }

    timeEndPeriod(1);

    getchar();
    return 0;
}