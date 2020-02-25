//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#ifndef PARTICLEEMITTER_H
#define PARTICLEEMITTER_H

#include "Math\Matrix.h"
#include "Math\Vect4D.h"
#include "Particle.h"

#include <list>

class ParticleEmitter
{
public:

	// big four
	ParticleEmitter();
	ParticleEmitter(const ParticleEmitter&) = default;
	~ParticleEmitter();
	ParticleEmitter& operator=(const ParticleEmitter&) = default;

	void update();
	void draw();
	void Execute(Vect4D& pos, Vect4D& vel, Vect4D& sc);

	//Matrix getFinalMatrix(const Matrix& scaleMatrix, const Matrix& transCamera, const Matrix& transParticle, const Matrix& rotParticle);
	//void matrixTransform(const Vect4D& scale, const Matrix& transCam, const Vect4D& transPar, const float& rot, Matrix& out);


private:
	Particle* pHead;
	Particle* buffer; // clean buffer to swap after each death

	Vect4D	vel_variance;
	Vect4D	pos_variance;
	Vect4D	start_position;
	Vect4D	start_velocity;

	float	max_life;
	int		max_particles;
	int		last_active_particle;
	float	spawn_frequency;
	float	last_spawn;
	float	last_loop;
	float	scale_variance;
};

#endif

// End of File
