//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Particle.h"



Particle::Particle()
{
	// construtor
	this->life = 0.0f;
	this->position.set(0.0f, 0.0f, 0.0f);
	this->velocity.set(0.0f, 0.0f, 0.0f);
	this->scale.set(1.0f, 1.0f, 1.0f);
	this->rotation = 0.0f;
	this->rotation_velocity = 0.25f;
}

Particle::Particle(const Particle& p)
{
	// copy the data only
	// this way of copying data is more efficient that element by element
	this->position = p.position;  // p.position
	this->velocity = p.velocity;
	this->scale = p.scale;
	this->rotation = p.rotation;
	this->rotation_velocity = p.rotation_velocity;
	this->life = p.life;
}

void Particle::CopyDataOnly(const Particle *p)
{
	// copy the data only
	// this way of copying data is more efficient that element by element
	this->position = p->position;
	this->velocity = p->velocity;
	this->scale = p->scale;
	this->rotation = p->rotation;
	this->rotation_velocity = p->rotation_velocity;
	this->life = p->life;
}

void Particle::Update(const float& time_elapsed)
{
	// amo: determinant always zero:
	//float MatrixScale = tmp.Determinant();  // !!! 4.8%

	// serious math below - magic secret sauce
	life += time_elapsed;
	

	//original ver position update
	position = position + (velocity * (time_elapsed));    //2.81%

	// !!! here
	//position.x += velocity.x * time_elapsed;
	//position.y += velocity.y * time_elapsed;
	//position.z += velocity.z * time_elapsed;

	//Vect4D tmp;
	//tmp.x = velocity.x * time_elapsed + position.x;
	//tmp.y = velocity.y * time_elapsed + position.y;
	//tmp.z = velocity.z* time_elapsed + position.z;

	//			x		y	  z
	// original ver.
	Vect4D axis(1.0f, 0.0f, 0.0f);
	Vect4D v(0.0f, 50.0f, 0.0f);
	position.Cross(axis, v);	//1.67%   

	//// cross ref
	////vout.x = (y * vin.z - z * vin.y);
	////vout.y = (z * vin.x - x * vin.z);
	////vout.z = (x * vin.y - y * vin.x);

	////v.x = (position.y * 0.0f - position.z * 0.0f);  = 0
	////v.y = (position.z * 1.0f- position.x * 0.0f);		= position.z
	////v.z = (position.x * 0.0f - position.y * 1.0f);	= position.y * -1.0f

	//v.x = 0.0f;
	//v.y = (position.z * 1.0f);
	//v.z = (position.y * -1.0f);

	// !!! here
	//float mag = sqrtf(position.z * position.z + position.y * position.y);
	//float magInv = 1.0f / mag;

	//float magTmp = sqrtf(tmp.z * tmp.z + tmp.y * tmp.y);
	//float magInvTmp = 1.0f / magTmp;



	// ref 
	//out.x = this->x * mag;
	//out.y = this->y * mag;
	//out.z = this->z * mag;

	// original
	v.norm(v);	//1.44%

	// cal ref
	//v.x = 0.0f;
	//v.y = (position.z * 1.0f) * mag;
	//v.y = position.z * mag;
	//v.z = (position.y * -1.0f) * mag;



	// original
	position = position + v * (0.07f * life);  // !!! 3.91%

	// ref
	//position.x = position.x + v.x * (0.07f * life);  // v.x is 0.0f stay same
	//position.y = position.y + v.x * (0.07f * life);

	// !!! here
	//position.y = position.y + position.z * (magInv * 0.07f * life);
	//position.z = position.z + position.y * (-1.0f * magInv * 0.07f * life);
	//position.y += position.z * (magInv * 0.07f * life);
	//position.z += position.y * (-1.0f * magInv * 0.07f * life);

	// tmp test
	//tmp.y = tmp.y + tmp.z * (magInvTmp * 0.07f * life);
	//tmp.z = tmp.z + tmp.y * (-1.0f * magInvTmp * 0.07f * life);

	//if (MatrixScale > 1.0f)
	//{
	//	MatrixScale = 1.0f / MatrixScale;
	//}

	rotation = rotation + rotation_velocity * time_elapsed * 4.0f;
	//rotation += rotation_velocity * time_elapsed * 4.0f;


	// after hard code position update
	//LoopTime : update : 0.638400 ms  draw : 83.068604 ms  tot : 83.707001
	//LoopTime : update : 0.716600 ms  draw : 83.915398 ms  tot : 84.631996
	//LoopTime : update : 0.731000 ms  draw : 87.958801 ms  tot : 88.689804
	//LoopTime : update : 0.632100 ms  draw : 86.919502 ms  tot : 87.551605
	//LoopTime : update : 0.652500 ms  draw : 82.473503 ms  tot : 83.125999
	//LoopTime : update : 0.651800 ms  draw : 84.108803 ms  tot : 84.760605
	//LoopTime : update : 0.787200 ms  draw : 86.181602 ms  tot : 86.968803
}


// End of file
