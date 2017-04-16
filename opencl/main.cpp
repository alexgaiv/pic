#include "mgl_suppress_warnings.h"
#include "common.h"
#include <stdlib.h>
#include <iostream>
#include <mgl2\mgl.h>
#include <time.h>
#include <grid.h>
#include "pic_kernel.h"
#include "tests.h"

#pragma comment(lib, "Winmm.lib")

using namespace cl;
using namespace std;

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename);

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

int main()
{
	timeBeginPeriod(1);
	
	//srand((unsigned)time(NULL));
	SetCurrentDirectory("C:\\Users\\HP\\Documents\\Visual Studio 2013\\Projects\\pic\\opencl");

	try
	{
		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_CPU);
		cout << "done\n";

		cl_Descriptor cld(ctx);

		YeeGrid grid0(Vector3f(0), Vector3f(1), Vector3i(16, 8, 8));
		YeeGrid grid1 = grid0;
		cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(16, 8, 8), Vector3i(4, 4, 4));
		//cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(1, 1, 1), Vector3i(1, 1, 1));
		PicKernel kernel(cld, "cl/kernel.cl", grid, true);

		Vector3f cellSize = grid.GetCellSize();
		Vector3i groupSize = grid.GetGroupSize();
		Vector3i numCells = grid.GetNumInnerCells();
		int totalNumCells = numCells.x * numCells.y * numCells.z;
		int groupNumCells = groupSize.x * groupSize.y * groupSize.z;

		cl_float3 *particles = new cl_float3[100 * totalNumCells];
		memset(particles, 0, 100 * totalNumCells * sizeof(cl_float3));

		FOR3(i, j, k, numCells)
		{
			int offset = 100 * idx(i, j, k, numCells);
			for (int p = 8; p < 100; p++)
			{
				float r = frand();
				particles[offset + p] = v2v(grid.GetMin() +
					(Vector3f(i, j, k) + Vector3f(r)) * cellSize);
			}

			particles[offset] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.1, 0.1, 0.1)) * cellSize);
			particles[offset + 1] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.1, 0.9, 0.1)) * cellSize);
			particles[offset + 2] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.1, 0.1, 0.9)) * cellSize);
			particles[offset + 3] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.1, 0.9, 0.9)) * cellSize);
			particles[offset + 4] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.9, 0.1, 0.1)) * cellSize);
			particles[offset + 5] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.9, 0.9, 0.1)) * cellSize);
			particles[offset + 6] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.9, 0.1, 0.9)) * cellSize);
			particles[offset + 7] = v2v(grid.GetMin() + (Vector3f(i, j, k) + Vector3f(0.9, 0.9, 0.9)) * cellSize);
		}

		int s1 = 12 * (groupSize.x + 2) * (groupSize.y + 2) * (groupSize.z + 2);
		int s2 = 12 * (numCells.x + 2) * (numCells.y + 2) * (numCells.z + 2);

		float *jx = new float[s2];
		float *jy = new float[s2];
		float *jz = new float[s2];
		cl_uint2 seed = {rand(), rand()};

		cl::Buffer jxMem(cld.ctx, CL_MEM_WRITE_ONLY, s2 * sizeof(float));
		cl::Buffer jyMem(cld.ctx, CL_MEM_WRITE_ONLY, s2 * sizeof(float));
		cl::Buffer jzMem(cld.ctx, CL_MEM_WRITE_ONLY, s2 * sizeof(float));
		cl::Buffer particlesMem(cld.ctx, CL_MEM_WRITE_ONLY, 100 * totalNumCells * sizeof(cl_float3));

		cld.queue.enqueueWriteBuffer(particlesMem, true, 0, 100 * totalNumCells * sizeof(cl_float3), particles);

		kernel.AddArg(s1 * sizeof(float), NULL);
		kernel.AddArg(s1 * sizeof(float), NULL);
		kernel.AddArg(s1 * sizeof(float), NULL);
		kernel.AddArg(jxMem);
		kernel.AddArg(jyMem);
		kernel.AddArg(jzMem);
		kernel.AddArg(particlesMem);
		kernel.AddArg(seed);

		kernel.Run();

		cld.queue.enqueueReadBuffer(jxMem, true, 0, s2 * sizeof(float), jx);
		cld.queue.enqueueReadBuffer(jyMem, true, 0, s2 * sizeof(float), jy);
		cld.queue.enqueueReadBuffer(jzMem, true, 0, s2 * sizeof(float), jz);
		grid.Jx.ReadBuffer();

		FOR3(i, j, k, grid.Jx.GetSize())
		{
			if (grid.Jx(i, j, k) == 6.0) __debugbreak();
		}

		Vector3i nc = numCells + Vector3i(2);
		FOR3(i, j, k, nc)
		{
			for (int p = 0; p < 12; p++) {
				if (jx[12 * idx(i, j, k, nc) + p] == 6.0f)
					__debugbreak();
			}

			if (i == 0 || j == 0 || k == 0 ||
				i == nc.x - 1 || j == nc.y - 1 || k == nc.z - 1)
			{
				for (int p = 0; p < 12; p++)
					jx[12 * idx(i, j, k, nc) + p] = 0.0f;
			}
		}

		Vector3i jxSize(3, 2, 2);
		Vector3i jySize(2, 3, 2);
		Vector3i jzSize(2, 2, 3);

		FOR3_(i, j, k, nc)
		{
			float j1 = 0.0f;
			float j2 = 0.0f;
			float j3 = 0.0f;
			float j4 = 0.0f;

			//if (i == 1 && j == 2 && k == 1) __debugbreak();
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

			grid0.Jx(i, j, k) = j1;
			grid0.Jx(i, j + 1, k) = j2;
			grid0.Jx(i, j, k + 1) = j3;
			grid0.Jx(i, j + 1, k + 1) = j4;
		}

		//FOR3(i, j, k, numCells)
		{
			//int offset = 100 * idx(i, j, k, numCells);
			int offset = 100 * idx(0, 1, 0, numCells);
			for (int t = 0; t < 1; t++)
			{
				cl_float4 coords = particles[offset + t];
				Particle p;
				p.coords = Vector3d(coords.x, coords.y, coords.z);
				p.momentum = Vector3d(0.0f);
				p.charge = electronCharge;
				p.mass = electronMass;
				p.factor = 1.0;
				grid1.DepositCurrents(p);
			}
		}

		FOR3_(i, j, k, grid0.Jx.GetSize())
		{
			float j1 = grid1.Jx(i, j, k);
			float j2 = grid.Jx(i, j, k);
			if (j1 == 0.0 && j2 == 0.0) continue;
			double err = abs(j1 - j2) / abs(j2);
			if (err > 0.0001) __debugbreak();
		}

		delete[] jx;
		delete[] jy;
		delete[] jz;
		delete[] particles;

		cout << "end";
	}
	catch (const Error &e)
	{
		cout << endl << e.what() << " failed, err=" << e.err();
	}

	timeEndPeriod(1);

	getchar();
	return 0;
}

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename)
{
	Vector3f vmin = grid.GetMin();
	Vector3f vmax = grid.GetMax();
	Vector3f cs = grid.GetCellSize();

	double *ex_plot = new double[latt.GetSize().x];
	double **ex_dens = new double*[latt.GetSize().y];
	for (int i = 0; i < latt.GetSize().y; i++)
		ex_dens[i] = new double[latt.GetSize().x];

	mglGraph gr;
	mglData y;

	Vector3d center = (vmin + vmax) / 2;

	for (int ix = 0; ix < latt.GetSize().x; ix++) {
		double x = vmin.x + (ix - 0.5) * cs.x;
		Vector3i s = latt.GetSize();
		ex_plot[ix] = latt(ix, s.y / 2, s.z / 2);
	}

	for (int ix = 0; ix < latt.GetSize().x; ix++) {
		for (int iy = 0; iy < latt.GetSize().y; iy++) {
			double x = vmin.x + (ix - 0.5) * cs.x;
			double y = vmin.y + (iy - 1.0) * cs.y;
			ex_dens[iy][ix] = latt(ix, iy, latt.GetSize().z / 2);
		}
	}

	const double A = 1111;
	gr.SetSize(1000, 401);
	gr.SubPlot(2, 1, 0);
	gr.SetOrigin(0, 0);
	gr.SetRanges(vmin.x, vmax.x, -A, A);
	gr.Label('x', "x");
	gr.Label('y', "E_x");
	gr.Axis();
	y.Set(ex_plot, latt.GetSize().x);
	gr.Plot(y);

	gr.SubPlot(2, 1, 1);
	gr.SetOrigin(0, 0);
	gr.SetRanges(vmin.x, vmax.x, vmin.y, vmax.y, -A, A);
	gr.Label('x', "x");
	gr.Label('y', "y");
	gr.Axis();
	gr.Box();
	gr.Colorbar();
	y.Set(ex_dens, latt.GetSize().y, latt.GetSize().x);
	gr.Dens(y);

	gr.WriteBMP(filename);

	delete[] ex_plot;
	for (int i = 0; i < latt.GetSize().y; i++)
		delete[] ex_dens[i];
	delete[] ex_dens;
}