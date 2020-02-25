//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#ifndef PARTICLE_H
#define PARTICLE_H

// include
#include "Vect4D.h"
#include "Matrix.h"
#include <xmmintrin.h>
#include <smmintrin.h>  

class Particle : Align16
{
public:
	friend class ParticleEmitter;

	// big four
	Particle();
	Particle(const Particle& copy);
	~Particle() = default;
	Particle &operator = (const Particle&) = default;

	void Update(const float& time_elapsed);
	void CopyDataOnly(const Particle *p); 


private:
	//Particle *next;
	//Particle *prev;

	Vect4D	position;
	Vect4D	velocity;
	Vect4D scale;
	//float	scale;

	float	life;
	float	rotation;
	float	rotation_velocity;
};


#endif 

// End of File
