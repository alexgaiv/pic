#ifndef _PARTICLE_MOVER_H_
#define _PARTICLE_MOVER_H_

#include "particle.h"

class ParticleMover
{
public:
	ParticleMover();
	ParticleMover(const Vector3d &E, const Vector3d &B)
		: E(E), B(B) { }

	Particle MoveParticle(const Particle &particle, double dt, int numSteps);
private:
	Vector3d E, B;

	double k;
	Vector3d u;
	double gamma;

	void iteration(Particle &pt, double dt);
};

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

Particle ParticleMover::MoveParticle(const Particle &particle, double dt, int numSteps)
{
	Particle pt = particle;
	const double mc = pt.mass * c;

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

#endif // _PARTICLE_MOVER_H_