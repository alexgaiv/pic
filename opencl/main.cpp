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

void BuildPlotEx(const cl_Lattice &latt, const cl_Grid &grid, const char *filename);

int main()
{
	timeBeginPeriod(1);

	try
	{
		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_GPU);
		cout << "done\n";

		cl_Descriptor cld(ctx);
		TestBoris(cld);
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