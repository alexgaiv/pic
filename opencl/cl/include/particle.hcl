#ifndef _PARTICLE_H_
#define _PARTICLE_H_

constant float c = 29979245800.0f;
constant float electronCharge = -4.80320427e-10f;
constant float electronMass = 9.10938215e-28f;

struct Particle
{
	float3 coords;
	float3 p; // momentum divided by (mass * c)
	float mass;
	float charge;
	float factor;

	float _invGamma; // 1 / sqrt(1 + p^2)
};

void particle_SetMomentum(struct Particle *pt, float3 momentum)
{
	pt->p = momentum / (pt->mass * c);
	pt->_invGamma = 1.0f / sqrt(1.0f + dot(pt->p, pt->p));
}

float3 particle_GetMomentum(struct Particle *pt) {
	return pt->p * (pt->mass * c);
}

void particle_SetVelocity(struct Particle *pt, float3 velocity)
{
	pt->p = velocity / sqrt(c * c - dot(velocity, velocity));
	pt->_invGamma = 1.0f / sqrt(1.0f + dot(pt->p, pt->p));
}

float3 particle_GetVelocity(struct Particle *pt)
{
	return pt->p * c * pt->_invGamma;
}

#endif // _PARTICLE_H_