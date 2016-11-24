#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "mgl_suppress_warnings.h"
#include <iostream>
#include <math.h>
#include <mgl2\mgl.h>
#include "utils.h"
#include "particle_mover.h"
#include "grid.h"
#include "tests.h"

using namespace std;

inline bool inrange(double a, double min, double max)
{
	return a >= min && a <= max;
}

void ColdPlasmaOscillations()
{
	YeeGrid grid(Vector3d(0), Vector3d(1, 0.125, 0.125), Vector3i(64, 8, 8));

	double A = 1111;
	double dt = 0.000000000000090453260945929121;
	int numSteps = 3688;
	double factor = 2941623.0;
	int dumpsPerIter = 16;

	Vector3d cs = grid.GetCellSize();
	Vector3d vmin = grid.GetMin();
	Vector3d vmax = grid.GetMax();
	double cellVolume = cs.x * cs.y * cs.z;

	Vector3i s = grid.Ex.GetSize();
	for (int i = 0; i < s.x; i++)
		for (int j = 0; j < s.y; j++)
			for (int k = 0; k < s.z; k++)
			{
				double x = vmin.x + (i - 0.5) * cs.x;
				grid.Ex(i, j, k) = A * cos(2 * M_PI * x);
			}

	vector<Particle> particles;

	s = grid.GetNumInnerCells();
	for (int i = 0; i < s.x; i++)
		for (int j = 0; j < s.y; j++)
			for (int k = 0; k < s.z; k++)
			{
				Vector3d lo = vmin + Vector3d(i, j, k) * cs;
				Vector3d hi = lo + cs;

				double cx = lo.x + cs.x*0.5;
				double f = 23133870163932 * (1.0 + 0.05 * sin(2 * M_PI * cx));
				double numParticles = f * cellVolume / factor;

				particles.reserve((unsigned)numParticles);

				for (int n = 0; n < numParticles; n++)
				{
					Particle p;
					p.charge = electronCharge;
					p.mass = electronMass;
					p.factor = factor;

					p.coords = Vector3d(
						frand(lo.x, hi.x),
						frand(lo.y, hi.y),
						frand(lo.z, hi.z));

					particles.push_back(p);
				}
			}

	double *ex_plot = new double[grid.Ex.GetSize().x];
	Vector3d center = (vmax + vmin) / 2;

	for (int i = 0; i < numSteps; i++)
	{
		cout << i << ' ';

		grid.ZeroiseJ();

		for (int k = 0, n = particles.size(); k < n; k++)
		{
			Particle &p = particles[k];
			if (inrange(p.coords.x, vmin.x, vmax.x) &&
				inrange(p.coords.y, vmin.y, vmax.y) &&
				inrange(p.coords.z, vmin.z, vmax.z))
			{
				grid.DepositCurrents(p);
			}
		}

		//grid.pbc_J();

		grid.SolveField(dt);

		for (int k = 0, n = particles.size(); k < n; k++)
		{
			Particle &p = particles[k];

			if (inrange(p.coords.x, vmin.x, vmax.x) &&
				inrange(p.coords.y, vmin.y, vmax.y) &&
				inrange(p.coords.z, vmin.z, vmax.z))
			{
				p = ParticleMover().MoveParticle(p, dt, 1, grid);
			}
		}

		if (i % dumpsPerIter == 0)
		{
			for (int x = 0; x < grid.Ex.GetSize().x; x++) {
				double x2 = vmin.x + (x - 0.5) * cs.x;
				ex_plot[x] = grid.InterpolateField(Vector3d(x2, center.y, center.z)).E.x;
			}
				
			mglGraph gr;
			mglData y;

			gr.Title("Ex");
			gr.SetOrigin(0, 0);
			gr.SetRanges(0.0, vmax.x, -A, A);
			gr.Label('x', "x");
			gr.Label('y', "Ex");
			gr.Axis();
			y.Set(ex_plot, grid.Ex.GetSize().x);
			gr.Plot(y);

			char filename[100];
			sprintf_s(filename, "plot/ColdPlasmaOscillations/ex-%d.bmp", i / dumpsPerIter + 1);
			gr.WriteBMP(filename);
		}
	}

	delete[] ex_plot;
}

int main()
{
	ColdPlasmaOscillations();
	//TestFdtd();

	getchar();
	return 0;
}