#include "tests.h"
#include "pic_kernel.h"
#include "mgl_suppress_warnings.h"
#include <mgl2\mgl.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace cl;

void BuildErrSubPlot(mglGraph *gr, const char *title, const char *legend, const char *color,
	float *data, int dataSize, float x_min, float x_max, float k_min, float k_max)
{
	mglData y;

	gr->Title(title);
	gr->ClearLegend();
	gr->AddLegend(legend, color);
	gr->SetRanges(x_min, x_max, data[dataSize - 1] * k_min, data[0] * k_max);
	gr->Label('x', "steps");
	gr->Label('y', "error");
	gr->Axis();
	gr->Legend();
	y.Set(data, dataSize);
	gr->Plot(y, color);
}

void BuildErrPlot(const char *filename, const char *title, float *data1, float *data2,
	int dataSize, float x_min, float x_max, float k_min, float k_max)
{
	mglGraph gr;

	gr.SetSize(1500, 600);
	gr.SubPlot(2, 1, 0);
	BuildErrSubPlot(&gr, title, "relative error", "b", data1, dataSize, x_min, x_max, k_min, k_max);
	gr.SubPlot(2, 1, 1);
	BuildErrSubPlot(&gr, title, "absolute error", "g", data2, dataSize, x_min, x_max, k_min, k_max);
	gr.WriteBMP(filename);
	gr.ClearFrame();
}

void TestGrid(cl_Descriptor &cld)
{
	const Vector3f vmin(0);
	const Vector3f vmax(1, 0.125, 0.125);
	const Vector3i numCells(64, 8, 8);
	const Vector3i groupNum(4, 4, 4);

	cl_Grid grid(cld, vmin, vmax, numCells, groupNum);
	const Vector3i groupSize = grid.GetGroupSize();

	int testSize = numCells.x * numCells.y * numCells.z;
	int *testData = new int[testSize];
	cl::Buffer testMem(cld.ctx, CL_MEM_WRITE_ONLY, testSize * sizeof(int));

	PicKernel kernel(cld, "cl/test_grid.cl", grid, true);
	kernel.AddArg(testMem);

	int dt = timeGetTime();
	kernel.Run();

	grid.Ex.ReadBuffer();
	grid.Ey.ReadBuffer();
	grid.Ez.ReadBuffer();
	grid.Bx.ReadBuffer();
	grid.By.ReadBuffer();
	grid.Bz.ReadBuffer();
	cld.queue.enqueueReadBuffer(testMem, true, 0, testSize * sizeof(int), testData);

	dt = timeGetTime() - dt;
	cout << "time: " << dt << "ms" << endl;

	bool passed = true;
	FOR3(i, j, k, numCells)
	{
		int idx = (k * numCells.y + j) * numCells.x + i;
		if (testData[idx] != 1) {
			passed = false;
			break;
		}
	}
	delete[] testData;
	cout << "grid test: " << (passed ? "passed" : "failed");
}

void TestBoris(cl_Descriptor &cld)
{
	cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(16, 8, 8), Vector3i(4, 4, 4));
	PicKernel kernel(cld, "cl/test_boris.cl", grid, true);

	const int startN = 100;
	const int dN = 50;
	const int numIters = 50;

	int maxN = dN * numIters;
	int memSize = numIters * sizeof(float);

	float r_rel_plot[numIters] = { };
	float r_abs_plot[numIters] = { };
	float p_rel_plot[numIters] = { };
	float p_abs_plot[numIters] = { };

	cl::Buffer buffers[8];

	kernel.AddArg(startN);
	kernel.AddArg(dN);
	kernel.AddArg(numIters);

	for (int i = 0; i < 8; i++) {
		buffers[i] = cl::Buffer(cld.ctx, CL_MEM_WRITE_ONLY, memSize);
		kernel.AddArg(buffers[i]);
	}

	kernel.Run();

	cld.queue.enqueueReadBuffer(buffers[0], true, 0, memSize, r_rel_plot);
	cld.queue.enqueueReadBuffer(buffers[1], true, 0, memSize, r_abs_plot);
	cld.queue.enqueueReadBuffer(buffers[2], true, 0, memSize, p_rel_plot);
	cld.queue.enqueueReadBuffer(buffers[3], true, 0, memSize, p_abs_plot);

	float p_abs[numIters] = { };
	ifstream file("1.bin", ios::binary);
	file.read((char *)p_abs, numIters * sizeof(float));

	int n = 0;
	for (int i = 0; i < numIters; i++) {
		if (p_abs[i] != p_abs_plot[i]) n++;
	}
	if (n != 0) __debugbreak();

	//ofstream fout("1.bin", ios::binary);
	//fout.write((char *)p_abs_plot, sizeof(float)*numIters);

	BuildErrPlot("plot/boris/test1-1.bmp", "Test 1 (Coords)",
		r_rel_plot, r_abs_plot, numIters, (float)startN, (float)maxN, 0.1f, 1.0f);
	BuildErrPlot("plot/boris/test1-2.bmp", "Test 1 (Momentum)",
		p_rel_plot, p_abs_plot, numIters, (float)startN, (float)maxN, 0.99f, 1.0f);

	cld.queue.enqueueReadBuffer(buffers[4], true, 0, memSize, r_rel_plot);
	cld.queue.enqueueReadBuffer(buffers[5], true, 0, memSize, r_abs_plot);
	cld.queue.enqueueReadBuffer(buffers[6], true, 0, memSize, p_rel_plot);
	cld.queue.enqueueReadBuffer(buffers[7], true, 0, memSize, p_abs_plot);

	BuildErrPlot("plot/boris/test2-1.bmp", "Test 2 (Coords)",
		r_rel_plot, r_abs_plot, numIters, (float)startN, (float)maxN, 1.0f, 1.0f);
	BuildErrPlot("plot/boris/test2-2.bmp", "Test 2 (Momentum)",
		p_rel_plot, p_abs_plot, numIters, (float)startN, (float)maxN, 1.0f, 1.0f);
}