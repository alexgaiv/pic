#ifndef _PARTICLE_H_
#define _PARTICLE_H_

constant float c = 29979245800.0f;
constant float electronCharge = -4.80320427e-10f;
constant float electronMass = 9.10938215e-28f;

struct Particle
{
	float3 coords;
	float3 momentum;
	float mass;
	float charge;
	float factor;
};

float3 particle_Velocity(struct Particle *p)
{
	float3 u = p->momentum / c;
	return p->momentum / sqrt(p->mass * p->mass + dot(u, u));
}

#endif // _PARTICLE_H_