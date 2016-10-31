#ifndef _PARTICLE_MOVER_H_
#define _PARTICLE_MOVER_H_

#include "particle.h"
#include "grid.h"

class ParticleMover
{
public:
	ParticleMover() {
		k = gamma = 0.0;
	}
	Particle MoveParticle(const Particle &particle, double dt, int numSteps, const YeeGrid &grid);
private:
	Vector3d E, B;
	double k;
	Vector3d u;
	double gamma;

	void iteration(Particle &pt, double dt);
};

#endif // _PARTICLE_MOVER_H_