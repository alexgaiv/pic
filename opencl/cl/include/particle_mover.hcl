#ifndef _PARTICLE_MOVER_H_
#define _PARTICLE_MOVER_H_

#include "particle.hcl"
#include "grid.hcl"

typedef struct
{
    float3 E, B;
    float k;
    float3 u;
    float gamma;
} ParticleMover;

void particleMover_Iteration(ParticleMover *pm, Particle *pt, float dt)
{
    float3 u_minus = pm->u + pm->k*pm->E;
    float3 t = (pm->k / pm->gamma) * pm->B;
    float3 s = 2 * t / (1 + dot(t, t));
    float3 u1 = u_minus + cross(u_minus, t);
    float3 u_plus = u_minus + cross(u1, s);

    pm->u = u_plus + pm->k*pm->E;
    pm->gamma = sqrt(1 + dot(pm->u, pm->u));

    float3 velocity = pm->u * c / pm->gamma;
    pt->coords += velocity * dt;
}

void MoveParticle(Particle *pt, Grid *grid, float dt)
{
    ParticleMover pm;
    
    float mc = pt->mass * c;
    FieldPoint f = grid_InterpolateField(grid, pt->coords);

    pm.E = f.E;
    pm.B = f.B;
    pm.k = pt->charge * dt / (2 * mc);
    pm.u = pt->p;
    pm.gamma = sqrt(1 + dot(pm.u, pm.u));

    particleMover_Iteration(&pm, pt, dt);
    pt->p = pm.u;
}

void MoveParticle2(Particle *pt, Grid *grid, float dt, int numSteps)
{
    ParticleMover pm;
    
    float mc = pt->mass * c;
    FieldPoint f = grid_InterpolateField(grid, pt->coords);

    pm.E = f.E;
    pm.B = f.B;
    pm.k = pt->charge * dt / (2 * mc);
    pm.u = pt->p;
    pm.gamma = sqrt(1 + dot(pm.u, pm.u));

    for (int i = 0; i < numSteps - 1; i++)
        particleMover_Iteration(&pm, pt, dt);

    dt *= 0.5;
    pm.k = pt->charge * dt / (2 * mc);

    particleMover_Iteration(&pm, pt, dt);
    float3 coords = pt->coords;
    particleMover_Iteration(&pm, pt, dt);

    pt->coords = coords;
    pt->p = pm.u;
}

#endif // _PARTICLE_MOVER_H_