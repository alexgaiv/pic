#include "mgl_suppress_warnings.h"
#include "common.h"
#include <stdlib.h>
#include <iostream>
#include <mgl2\mgl.h>
#include "pic_kernel.h"
#include "tests.h"

#pragma comment(lib, "Winmm.lib")

using namespace cl;
using namespace std;

ostream &operator<<(ostream &os, const cl_float3 &v)
{
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename);

void TestBoris(cl_Descriptor &cld)
{
	cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(16, 16, 16), Vector3i(4, 4, 4));
	PicKernel kernel(cld, "cl/kernel.cl", grid, true);

	cl_float2 err[2] = { };
	cl_float3 rp[2];

	cl::Buffer errMem(cld.ctx, CL_MEM_WRITE_ONLY, sizeof(err));
	cl::Buffer rpMem(cld.ctx, CL_MEM_WRITE_ONLY, sizeof(rp));
	kernel.AddArg(errMem);
	kernel.AddArg(rpMem);

	kernel.Run();

	cld.queue.enqueueReadBuffer(errMem, true, 0, sizeof(err), err);
	cld.queue.enqueueReadBuffer(rpMem, true, 0, sizeof(rp), rp);

	cout << "r = " << rp[0] << ", p = " << rp[1] << endl;

	cout << "r_error_abs = " << err[0].x << ", p_error_abs = " << err[1].x << endl;
	cout << "r_error_rel = " << err[0].y << ", p_error_rel = " << err[1].y << "\n\n";
}

int main()
{
	timeBeginPeriod(1);

	try
	{
		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_CPU);
		cout << "done\n";

		cl_Descriptor cld(ctx);
		TestBoris(cld);
	}
	catch (const Error &e)
	{
		cout << endl << e.what() << " failed, err=" << e.err();
	}

	//BuildPlotEx(grid.Ex, grid, "plot/ex.bmp");

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