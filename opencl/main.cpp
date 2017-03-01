#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#define __CL_ENABLE_EXCEPTIONS
#define NOMINMAX

#include "mgl_suppress_warnings.h"
#include <Windows.h>
#include <stdlib.h>
#include <iostream>
#include <cl\cl.hpp>
#include <mgl2\mgl.h>
#include "mathtypes.h"
#include "utils.h"
#include "cl_grid.h"

#pragma comment(lib, "Winmm.lib")

using namespace cl;
using namespace std;

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename);

int main()
{
	timeBeginPeriod(1);

	try
	{
		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_CPU);
		cout << "done\n";

		Device device = ctx.getInfo<CL_CONTEXT_DEVICES>()[0];
		CommandQueue queue(ctx, device);
		Program program = CreateProgramFromFile("kernel.cl", ctx);

		try
		{
			cout << "building a program...";
			program.build("-I cl");
			cout << "done\n";
			string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
			if (log != "") OutputDebugStringA(log.c_str());
		}
		catch (const Error &e)
		{
			string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
			OutputDebugStringA(log.c_str());
			throw e;
		}

		cl_Grid grid(queue, Vector3f(0), Vector3f(1, 0.125, 0.125), Vector3i(64, 8, 8));
		Vector3i groupNum(4, 4, 4);
		Vector3i numCells = grid.GetNumInnerCells();
		Vector3i groupSize = numCells / groupNum;

		Kernel kernel(program, "main");
		kernel.setArg(0, v2v(grid.GetMin()));
		kernel.setArg(1, v2v(grid.GetMax()));
		kernel.setArg(2, v2v(grid.GetNumInnerCells()));

		kernel.setArg(3, grid.Ex.buffer);
		kernel.setArg(4, grid.Ey.buffer);
		kernel.setArg(5, grid.Ez.buffer);
		kernel.setArg(6, grid.Bx.buffer);
		kernel.setArg(7, grid.By.buffer);
		kernel.setArg(8, grid.Bz.buffer);
		kernel.setArg(9, grid.Jx.buffer);
		kernel.setArg(10, grid.Jy.buffer);
		kernel.setArg(11, grid.Jz.buffer);

		kernel.setArg(12, groupNum.x * (groupNum.y + 1) * (groupNum.z + 1) * sizeof(real_t), NULL);
		kernel.setArg(13, (groupNum.x + 1) * groupNum.y * (groupNum.z + 1) * sizeof(real_t), NULL);
		kernel.setArg(14, (groupNum.x + 1) * (groupNum.y + 1) * groupNum.z * sizeof(real_t), NULL);

		kernel.setArg(15, (groupNum.x + 1) * groupNum.y * groupNum.z * sizeof(real_t), NULL);
		kernel.setArg(16, groupNum.x * (groupNum.y + 1) * groupNum.z * sizeof(real_t), NULL);
		kernel.setArg(17, groupNum.x * groupNum.y * (groupNum.z + 1) * sizeof(real_t), NULL);

		kernel.setArg(18, groupNum.x * (groupNum.y + 1) * (groupNum.z + 1) * sizeof(real_t), NULL);
		kernel.setArg(19, (groupNum.x + 1) * groupNum.y * (groupNum.z + 1) * sizeof(real_t), NULL);
		kernel.setArg(20, (groupNum.x + 1) * (groupNum.y + 1) * groupNum.z * sizeof(real_t), NULL);

		int dt = timeGetTime();
		queue.enqueueNDRangeKernel(kernel, NullRange,
			NDRange(numCells.x, numCells.y, numCells.z),
			NDRange(groupSize.x, groupSize.y, groupSize.z));

		grid.Ex.ReadBuffer();
		grid.Ey.ReadBuffer();
		grid.Ez.ReadBuffer();
		grid.Bx.ReadBuffer();
		grid.By.ReadBuffer();
		grid.Bz.ReadBuffer();

		dt = timeGetTime() - dt;
		cout << dt << "ms" << endl;

		BuildPlotEx(grid.Ex, grid, "plot/ex.bmp");
		
		queue.flush();
		queue.finish();
	}
	catch (const Error &e)
	{
		cout << endl << e.what() << " failed, err=" << e.err();
	}
	
	cout << "\nend";
	getchar();
	timeEndPeriod(1);
	return 0;
}

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename)
{
	Vector3f vmin = grid.GetMin();
	Vector3f vmax = grid.GetMax();
	Vector3f cs = grid.GetCellSize();

	double *ex_plot = new double[latt.GetSize().x - 2];
	double **ex_dens = new double*[latt.GetSize().x - 2];
	for (int i = 0; i < latt.GetSize().x - 2; i++)
		ex_dens[i] = new double[latt.GetSize().y - 2];

	mglGraph gr;
	mglData y;

	Vector3d center = (vmin + vmax) / 2;

	for (int ix = 1; ix < latt.GetSize().x - 1; ix++) {
		double x = vmin.x + (ix - 0.5) * cs.x;
		Vector3i s = latt.GetSize();
		ex_plot[ix - 1] = latt(ix, s.y / 2, s.z / 2);
	}

	for (int ix = 1; ix < latt.GetSize().x - 1; ix++) {
		for (int iy = 1; iy < latt.GetSize().y - 1; iy++) {
			double x = vmin.x + (ix - 0.5) * cs.x;
			double y = vmin.y + (iy - 1.0) * cs.y;
			ex_dens[ix - 1][iy - 1] = latt(ix, iy, latt.GetSize().z / 2);
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
	y.Set(ex_plot, latt.GetSize().x - 2);
	gr.Plot(y);

	gr.SubPlot(2, 1, 1);
	gr.SetOrigin(0, 0);
	gr.SetRanges(vmin.x, vmax.x, vmin.y, vmax.y, -A, A);
	gr.Label('x', "x");
	gr.Label('y', "y");
	gr.Axis();
	gr.Box();
	gr.Colorbar();
	y.Set(ex_dens, latt.GetSize().x - 2, latt.GetSize().y - 2);
	gr.Dens(y);

	gr.WriteBMP(filename);

	delete[] ex_plot;
	for (int i = 0; i < latt.GetSize().x - 2; i++)
		delete[] ex_dens[i];
	delete[] ex_dens;
}