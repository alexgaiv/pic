#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <math.h>
#include <mgl2\mgl.h>
#include "particle_mover.h"
#include "grid.h"

using namespace std;

struct ErrorStruct
{
	double abs;
	double rel;
};

bool CmpReal(double a, double b) {
	return abs(a - b) <= 0.0001;
}

ostream &operator<<(ostream &os, const Vector3d &v)
{
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}

void ConstantField(YeeGrid &grid, const Vector3d &E, const Vector3d &B)
{
	fill(grid.Ex.GetData().begin(), grid.Ex.GetData().end(), E.x);
	fill(grid.Ey.GetData().begin(), grid.Ey.GetData().end(), E.y);
	fill(grid.Ez.GetData().begin(), grid.Ez.GetData().end(), E.z);

	fill(grid.Bx.GetData().begin(), grid.Bx.GetData().end(), B.x);
	fill(grid.By.GetData().begin(), grid.By.GetData().end(), B.y);
	fill(grid.Bz.GetData().begin(), grid.Bz.GetData().end(), B.z);
}

void TestBoris_1(int steps, ErrorStruct &r_error, ErrorStruct &p_error, bool output = false)
{
	Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;

	const double mc = p.mass * c;
	const double E0 = 43;

	int n1 = output ? 3 : 1;
	int n2 = output ? -1 : 0;
	for (int i = 0; i < n1; i++) // for each axis
	{
		for (int k = 1; k >= n2; k -= 2) // for each direction
		{
			Vector3d E;
			Vector3d B;

			double Ek = E0 * k;
			E[i] = Ek;

			YeeGrid grid(Vector3d(-30), Vector3d(30), Vector3i(15, 20, 12));
			ConstantField(grid, E, B);

			double dt = mc / (p.charge * Ek * steps);

			Particle p1 = ParticleMover().MoveParticle(p, dt, steps, grid);

			Vector3d r_theor;
			r_theor[i] = mc * c / (p.charge * Ek) * (sqrt(2.0) - 1.0);
			Vector3d p_theor;
			p_theor[i] = mc;

			r_error.abs = (p1.coords - r_theor).Length();
			p_error.abs = (p1.momentum - p_theor).Length();
			r_error.rel = r_error.abs / r_theor.Length();
			p_error.rel = p_error.abs / p_theor.Length();

			if (output)
			{
				cout << "Test 1 (N = " << steps << ")\n";
				cout << "E = " << E << endl;
				cout << "r = " << r_theor << ", p = " << p_theor << endl;
				cout << "r = " << p1.coords << ", p = " << p1.momentum << endl;
				cout << "r_error_abs = " << r_error.abs << ", p_error_abs = " << p_error.abs << endl;
				cout << "r_error_rel = " << r_error.rel << ", p_error_rel = " << p_error.rel << "\n\n";
			}
		}
	}
}

void TestBoris_2(int steps, ErrorStruct &r_error, ErrorStruct &p_error, bool output = false)
{
	Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;

	const double mc = p.mass * c;
	const double p0 = 5;
	const double B0 = 57;

	Vector3d E;
	Vector3d B(0, 0, B0);
	p.momentum = Vector3d(p0, 0, 0);

	YeeGrid grid(Vector3d(-30), Vector3d(30), Vector3i(15, 20, 12));
	ConstantField(grid, E, B);

	double dt = M_PI * mc / (abs(p.charge) * B0 * steps) * sqrt(1 + pow(p0 / mc, 2));

	Particle p1 = ParticleMover().MoveParticle(p, dt, steps, grid);

	Vector3d r_theor(0, -2 * p0*c / (p.charge * B0), 0);
	Vector3d p_theor(-p0, 0, 0);

	r_error.abs = (p1.coords - r_theor).Length();
	p_error.abs = (p1.momentum - p_theor).Length();
	r_error.rel = r_error.abs / r_theor.Length();
	p_error.rel = p_error.abs / p_theor.Length();

	if (output)
	{
		cout << "\n\nTest 2 (N = " << steps << ")\n";
		cout << "B = " << B << endl;
		cout << "r = " << r_theor << ", p = " << p_theor << endl;
		cout << "r = " << p1.coords << ", p = " << p1.momentum << endl;
		cout << "r_error_abs = " << r_error.abs << ", p_error_abs = " << p_error.abs << endl;
		cout << "r_error_rel = " << r_error.rel << ", p_error_rel = " << p_error.rel << "\n\n";
	}
}

void BuildPlots()
{
	ErrorStruct r_error = { };
	ErrorStruct p_error = { };

	mglData y;
	mglGraph gr;

	int N;
	const int startN = 100;
	const int dN = 50;
	const int maxN = 2500;
	const int s = maxN / dN;

	double r_rel_plot[s] = { };
	double r_abs_plot[s] = { };
	double p_rel_plot[s] = { };
	double p_abs_plot[s] = { };

	/*Test 1*/
	N = startN;
	for (int i = 0; i < s; i++, N += dN)
	{
		TestBoris_1(N, r_error, p_error);
		r_rel_plot[i] = r_error.rel;
		r_abs_plot[i] = r_error.abs;
		p_rel_plot[i] = p_error.rel;
		p_abs_plot[i] = p_error.abs;
	}

	gr.Title("Test 1 (Coords)");
	gr.AddLegend("relative error", "b");
	gr.AddLegend("absolute error", "g");
	gr.SetRanges(startN, maxN, 0, r_rel_plot[0] * 1.1);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	gr.Legend();
	y.Set(r_rel_plot, s);
	gr.Plot(y, "b");
	y.Set(r_abs_plot, s);
	gr.Plot(y, "g");

	gr.WriteBMP("plot/test1-1.bmp");
	gr.ClearFrame();

	gr.SetSize(1500, 600);
	gr.SubPlot(2, 1, 0);
	gr.Title("Test 1 (Momentum, rel)");
	gr.SetRanges(startN, maxN, 0, 1e-13);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	y.Set(p_rel_plot, s);
	gr.Plot(y);

	gr.SubPlot(2, 1, 1);
	gr.Title("Test 1 (Momentum, abs)");
	gr.SetRanges(startN, maxN, 0, 1e-29);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	y.Set(p_abs_plot, s);
	gr.Plot(y, "g");

	gr.WriteBMP("plot/test1-2.bmp");
	gr.ClearFrame();

	/*Test 2*/
	N = startN;
	for (int i = 0; i < s; i++, N += dN)
	{
		TestBoris_2(N, r_error, p_error);
		r_rel_plot[i] = r_error.rel;
		r_abs_plot[i] = r_error.abs;
		p_rel_plot[i] = p_error.rel;
		p_abs_plot[i] = p_error.abs;
	}

	gr.SubPlot(2, 1, 0);
	gr.Title("Test 2 (Coords, rel)");
	gr.SetRanges(startN, maxN, 0, r_rel_plot[0] * 1.1);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	y.Set(r_rel_plot, s);
	gr.Plot(y);

	gr.SubPlot(2, 1, 1);
	gr.Title("Test 2 (Coords, abs)");
	gr.SetRanges(startN, maxN, 0, r_abs_plot[0] * 1.1);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	y.Set(r_abs_plot, s);
	gr.Plot(y, "g");

	gr.WriteBMP("plot/test2-1.bmp");
	gr.ClearFrame();

	gr.SetSize(600, 400);
	gr.SubPlot(1, 1, 0);
	gr.Title("Test 2 (Momentum)");
	gr.AddLegend("relative error", "b");
	gr.AddLegend("absolute error", "g");
	gr.SetRanges(startN, maxN, 0, 1e-4);
	gr.Label('x', "steps");
	gr.Label('y', "error");
	gr.Axis();
	gr.Legend();
	y.Set(p_rel_plot, s);
	gr.Plot(y, "b");
	y.Set(p_abs_plot, s);
	gr.Plot(y, "g");

	gr.WriteBMP("plot/test2-2.bmp");
}

void TestBoris()
{
	ErrorStruct r_error = { };
	ErrorStruct p_error = { };

	TestBoris_1(100, r_error, p_error, true);
	TestBoris_2(100, r_error, p_error, true);

	cout << "building plots...";
	BuildPlots();
	cout << "done.";
}

void TestGrid()
{
	YeeGrid grid(Vector3d(-10, -20, 15), Vector3d(30, 20, 25), Vector3i(15, 10, 20));

	Vector3d cs = grid.GetCellSize();
	Vector3d vmin = grid.GetMin();
	Vector3d vmax = grid.GetMax();
	Vector3d step = (vmax - vmin) / 10;

	bool passed = true;
	Vector3i s = grid.Ex.GetSize();
	for (int i = 0; i < s.x; i++)
		for (int j = 0; j < s.y; j++)
			for (int k = 0; k < s.z; k++)
				grid.Ex(i, j, k) = vmin.x + (i - 0.5) * cs.x;

	for (double x = vmin.x; x <= vmax.x; x += step.x)
		for (double y = vmin.y; y <= vmax.y; y += step.y)
			for (double z = vmin.z; z <= vmax.z; z += step.z)
			{
				double Ex = grid.InterpolateField(Vector3d(x, y, z)).E.x;
				if (!CmpReal(x, Ex)) {
					passed = false;
					goto exit_loop1;
				}
			}
exit_loop1:
	passed &= CmpReal(vmax.x, grid.InterpolateField(vmax).E.x);
	cout << "grid test (Ex): " << (passed ? "passed" : "failed") << endl;


	passed = true;
	s = grid.Bx.GetSize();
	for (int i = 0; i < s.x; i++)
		for (int j = 0; j < s.y; j++)
			for (int k = 0; k < s.z; k++)
				grid.Bx(i, j, k) = vmin.x + i*cs.x;

	for (double x = vmin.x; x <= vmax.x; x += step.x)
		for (double y = vmin.y; y <= vmax.y; y += step.y)
			for (double z = vmin.z; z <= vmax.z; z += step.z)
			{
				double Bx = grid.InterpolateField(Vector3d(x, y, z)).B.x;
				if (!CmpReal(x, Bx)) {
					passed = false;
					goto exit_loop2;
				}
			}
exit_loop2:
	passed &= CmpReal(vmax.x, grid.InterpolateField(vmax).B.x);
	cout << "grid test (Bx): " << (passed ? "passed" : "failed") << endl;
}

int main()
{
	TestGrid();
	getchar();
	return 0;
}