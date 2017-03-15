#ifndef _PARTICLE_H_
#define _PARTICLE_H_

constant float c = 29979245800.0;
constant float electronCharge = -4.80320427e-10;
constant float electronMass = 9.10938215e-28;

struct Particle
{
	float3 coords;
	float3 momentum;
	float mass;
	float charge;
	float factor;
};

#endif // _PARTICLE_H_