#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "mathtypes.h"

const double c = 29979245800.0;
const double electronCharge = -4.80320427e-10;
const double electronMass = 9.10938215e-28;

class Particle
{
public:
	Vector3d coords;
	Vector3d momentum;
	double mass;
	double charge;

	Particle() {
		mass = charge = 0.0;
	}

	Vector3d Velocity() const {
		return momentum / sqrt(mass * mass + (momentum / c).Square());
	}
};

#endif