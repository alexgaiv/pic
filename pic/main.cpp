#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "mgl_suppress_warnings.h"
#include <iostream>
#include <math.h>
#include <mgl2\mgl.h>
#include "utils.h"
#include "array3D.h"
#include "particle_mover.h"
#include "grid.h"
#include "ensemble.h"
#include "tests.h"

using namespace std;

void ColdPlasmaOscillations()
{
	YeeGrid grid(Vector3d(0), Vector3d(1, 0.125, 0.125), Vector3i(64, 8, 8));

	const double A = 1111;
	const double dt = 0.000000000000090453260945929121;
	const int numSteps = 3688;
	const double factor = 2941623;
	const int dumpsPerIter = 16;

	Vector3i numCells = grid.GetNumInnerCells();
	Vector3d cs = grid.GetCellSize();
	Vector3d vmin = grid.GetMin();
	Vector3d vmax = grid.GetMax();
	Vector3d size = vmax - vmin;
	Vector3d center = (vmin + vmax) / 2;
	double cellVolume = cs.x * cs.y * cs.z;

	const Vector3i es(3, 3, 3);
	Array3D<Ensemble> ensembles(es);
	
	FOR3(ei, ej, ek, es)
	{
		Vector3i s = numCells / es;
		Vector3d min = vmin + Vector3d(ei, ej, ek) * Vector3d(s) * cs;

		if (ei == es.x - 1) s.x += numCells.x % es.x;
		if (ej == es.y - 1) s.y += numCells.y % es.y;
		if (ek == es.z - 1) s.z += numCells.z % es.z;
		
		ensembles(ei, ej, ek) = Ensemble(min, cs, s);
	}

	double *ex_plot = new double[grid.Ex.GetSize().x];
	double **ex_dens = new double*[grid.Ex.GetSize().y];
	for (int i = 0; i < grid.Ex.GetSize().y; i++)
		ex_dens[i] = new double[grid.Ex.GetSize().x];

	Vector3i s = grid.Ex.GetSize();
	FOR3(i, j, k, s)
	{
		double x = vmin.x + (i - 0.5) * cs.x;
		grid.Ex(i, j, k) = A * cos(2 * M_PI * x);
	}

	FOR3(ei, ej, ek, es)
	{
		Ensemble &e = ensembles(ei, ej, ek);
		Vector3i s = e.GetSize();

		FOR3(i, j, k, s)
		{
			Vector3d lo = e.GetMin() + Vector3d(i, j, k) * cs;
			Vector3d hi = lo + cs;

			double cx = lo.x + cs.x*0.5;
			double f = 23133870163932 * (1.0 + 0.05 * sin(2 * M_PI * cx));
			int numParticles = int(f * cellVolume / factor);

			ParticleSystem &system = (e(i, j, k) = ParticleSystem(numParticles));

			for (int q = 0 ; q < numParticles; q++)
			{
				Particle p;
				p.charge = electronCharge;
				p.mass = electronMass;
				p.factor = factor;

				p.coords = Vector3d(
					frand(lo.x, hi.x),
					frand(lo.y, hi.y),
					frand(lo.z, hi.z));

				system.Add(p);
			}
		}
	}
	
	for (int step = 0; step < numSteps; step++)
	{
		cout << step << ' ';

		grid.ZeroiseJ();

		FOR3(ei, ej, ek, es)
		{
			Ensemble &e = ensembles(ei, ej, ek);
			Vector3i s = e.GetSize();
			FOR3(i, j, k, s)
			{
				ParticleSystem &system = e(i, j, k);
				system.Begin();
				while (!system.End())
				{
					grid.DepositCurrents(system.Current());
					system.Next();
				}
			}
		}

		grid.PbcJ();

		grid.SolveField(dt);

		FOR3(ei, ej, ek, es)
		{
			Ensemble &e = ensembles(ei, ej, ek);
			Vector3i s = e.GetSize();
			FOR3(i, j, k, s)
			{
				ParticleSystem &system = e(i, j, k);
				system.Begin();
				while (!system.End())
				{
					Particle &p = system.Current();
					p = ParticleMover().MoveParticle(p, dt, 1, grid);

					if (p.coords.x < vmin.x) p.coords.x += size.x;
					else if (p.coords.x > vmax.x) p.coords.x -= size.x;
					if (p.coords.y < vmin.y) p.coords.y += size.y;
					else if (p.coords.y > vmax.y) p.coords.y -= size.y;
					if (p.coords.z < vmin.z) p.coords.z += size.z;
					else if (p.coords.z > vmax.z) p.coords.z -= size.z;

					system.Next();
				}
			}
		}

		if (step % dumpsPerIter == 0)
		{
			mglGraph gr;
			mglData y;

			for (int ix = 0; ix < grid.Ex.GetSize().x; ix++) {
				double x = vmin.x + (ix - 0.5) * cs.x;
				ex_plot[ix] = grid.InterpolateField(Vector3d(x, center.y, center.z)).E.x;
			}

			for (int ix = 0; ix < grid.Ex.GetSize().x; ix++) {
				for (int iy = 0; iy < grid.Ex.GetSize().y; iy++) {
					double x = vmin.x + (ix - 0.5) * cs.x;
					double y = vmin.y + (iy - 1.0) * cs.y;
					ex_dens[iy][ix] = grid.InterpolateField(Vector3d(x, y, center.z)).E.x;
				}
			}
			
			gr.SetSize(1000, 401);
			gr.SubPlot(2, 1, 0);
			gr.SetOrigin(0, 0);
			gr.SetRanges(vmin.x, vmax.x, -A, A);
			gr.Label('x', "x");
			gr.Label('y', "E_x");
			gr.Axis();
			y.Set(ex_plot, grid.Ex.GetSize().x);
			gr.Plot(y);

			gr.SubPlot(2, 1, 1);
			gr.SetOrigin(0, 0);
			gr.SetRanges(vmin.x, vmax.x, vmin.y, vmax.y, -A, A);
			gr.Label('x', "x");
			gr.Label('y', "y");
			gr.Axis();
			gr.Box();
			gr.Colorbar();
			y.Set(ex_dens, grid.Ex.GetSize().y, grid.Ex.GetSize().x);
			gr.Dens(y);

			char filename[100];
			sprintf_s(filename, "plot/ColdPlasmaOscillations/ex-%d.bmp", step / dumpsPerIter + 1);
			gr.WriteBMP(filename);
		}
	}

	delete[] ex_plot;
	for (int i = 0; i < grid.Ex.GetSize().y; i++)
		delete[] ex_dens[i];
	delete[] ex_dens;
}

int main()
{
	//ColdPlasmaOscillations();
	//TestBoris();
	TestGrid();
	getchar();
	return 0;
}