#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <math.h>
#include "particle_mover.h"

using namespace std;

ostream &operator<<(ostream &os, const Vector3d &v)
{
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}

void Test1()
{
	Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;

	const int N = 100;
	const double mc = p.mass * c;
	const double E0 = 43;

	cout << "\n\nTest 1 (N = " << N << ")\n\n";

	for (int i = 0; i < 3; i++) // for each axis
	{
		for (int k = 1; k >= -1; k -= 2) // for each direction
		{
			Vector3d E;
			Vector3d B;

			double Ek = E0 * k;
			E[i] = Ek;

			double dt = mc / (p.charge * Ek * N);

			Particle p1 = ParticleMover(E, B).MoveParticle(p, dt, N);

			Vector3d r_theor;
			r_theor[i] = mc * c / (p.charge * Ek) * (sqrt(2.0) - 1.0);
			Vector3d p_theor;
			p_theor[i] = mc;

			double r_error = (p1.coords - r_theor).Length();
			double p_error = (p1.momentum - p_theor).Length();

			cout << "E = " << E << endl;
			cout << "r = " << r_theor << ", p = " << p_theor << endl;
			cout << "r = " << p1.coords << ", p = " << p1.momentum << endl;
			cout << "r_error = " << r_error << ", p_error = " << p_error << "\n\n";
		}
	}
}

void Test2()
{
	Particle p;
	p.mass = electronMass;
	p.charge = electronCharge;

	const int N = 100;
	const double mc = p.mass * c;
	const double p0 = 5;
	const double B0 = 57;

	Vector3d E;
	Vector3d B(0, 0, B0);
	p.momentum = Vector3d(p0, 0, 0);

	double dt = M_PI * mc / (abs(p.charge) * B0 * N) * sqrt(1 + pow(p0 / mc, 2));

	Particle p1 = ParticleMover(E, B).MoveParticle(p, dt, N);

	Vector3d r_theor(0, -2 * p0*c / (p.charge * B0), 0);
	Vector3d p_theor(-p0, 0, 0);

	double r_error = (p1.coords - r_theor).Length();
	double p_error = (p1.momentum - p_theor).Length();

	cout << "\n\nTest 2 (N = " << N << ")\n\n";
	cout << "B = " << B << endl;
	cout << "r = " << r_theor << ", p = " << p_theor << endl;
	cout << "r = " << p1.coords << ", p = " << p1.momentum << endl;
	cout << "r_error = " << r_error << ", p_error = " << p_error << "\n\n";
}

int main()
{
	Test1();
	Test2();
	getchar();
	return 0;
}