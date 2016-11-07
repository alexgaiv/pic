#include "particle_mover.h"

void ParticleMover::iteration(Particle &pt, double dt)
{
	Vector3d u_minus = u + k*E;
	Vector3d t = (k / gamma) * B;
	Vector3d s = 2.0 * t / (1 + t.Square());
	Vector3d u1 = u_minus + Cross(u_minus, t);
	Vector3d u_plus = u_minus + Cross(u1, s);

	u = u_plus + k*E;
	gamma = sqrt(1 + u.Square());

	Vector3d velocity = u * c / gamma;
	pt.coords += velocity * dt;
}

Particle ParticleMover::MoveParticle(const Particle &particle, double dt, int numSteps, const YeeGrid &grid)
{
	Particle pt = particle;
	const double mc = pt.mass * c;

	FieldPoint f = grid.InterpolateField(pt.coords);
	E = f.E;
	B = f.B;

	k = pt.charge * dt / (2 * mc);
	u = pt.momentum / mc;
	gamma = sqrt(1 + u.Square());

	for (int i = 0; i < numSteps - 1; i++)
		iteration(pt, dt);

	dt *= 0.5;
	k = pt.charge * dt / (2 * mc);

	iteration(pt, dt);
	Vector3d coords = pt.coords;
	iteration(pt, dt);

	pt.coords = coords;
	pt.momentum = mc * u;
	return pt;
}