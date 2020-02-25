//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "DO_NOT_MODIFY\OpenGLInterface.h"
#include "ParticleEmitter.h"
#include "Settings.h"


extern PerformanceTimer GlobalTimer;
extern Matrix inverseCameraMatrix;

static const unsigned char squareColors[] =
{
	// ----------------------------
	//  point is actually a quad   
	//  set color on each vertex   
	// ----------------------------
	//    Vert A = Yellow 
	//    Vert B = Yellow
	//    Vert C = Yellow
	//    Vert D = Yellow
	// ----------------------------

	200,  200,  0,  255,
	200,  200,  0,  255,
	200,  200,  0,  255,
	200,  200,  0,  255,
};

static const float squareVertices[] =
{
	// --------------------
	//   Vert A = (x,y,z)
	//   Vert B = (x,y,z)
	//   Vert C = (x,y,z)
	//   Vert D = (x,y,z)
	// --------------------

	-0.015f,  -0.015f, 0.0f, // Size of Triangle
	 0.015f,  -0.015f, 0.0f, // Size of Triangle
	-0.015f,   0.015f, 0.0f, // Size of Triangle
	 0.015f,   0.015f, 0.0f, // Size of Triangle
};

ParticleEmitter::ParticleEmitter()
	: start_position(0.0f, 0.0f, 0.0f),
	start_velocity(0.0f, 1.0f, 0.0f),
	max_life(MAX_LIFE),
	max_particles(NUM_PARTICLES),
	spawn_frequency(0.0000001f),
	last_active_particle(-1),
	vel_variance(1.0f, 4.0f, 0.4f),
	pos_variance(1.0f, 1.0f, 1.0f),
	scale_variance(2.5f)
{
	// OpenGL goo... don't worrry about this
	glVertexPointer(3, GL_FLOAT, 0, squareVertices);
	glEnableClientState(GL_VERTEX_ARRAY);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, squareColors);
	glEnableClientState(GL_COLOR_ARRAY);

	pHead = new Particle[NUM_PARTICLES];
	Particle* t = pHead;
	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		t->life = 0.0f;
		t->position = start_position;
		t->velocity = start_velocity;
		t->scale.set(1.0f, 1.0f, 1.0f, 1.0f);
		this->Execute(t->position, t->velocity, t->scale);
		t++;
	}

	buffer = new Particle[NUM_PARTICLES];
	memcpy(buffer, pHead, NUM_PARTICLES * sizeof(Particle));
	
	GlobalTimer.Toc();

	last_spawn = GlobalTimer.TimeInSeconds();
	last_loop = GlobalTimer.TimeInSeconds();
}

ParticleEmitter::~ParticleEmitter()
{
	delete[] pHead;
	delete[] buffer;
}

void ParticleEmitter::update()
{
	// get current time
	GlobalTimer.Toc();
	float current_time = GlobalTimer.TimeInSeconds();

	// spawn particles
	float time_elapsed = current_time - last_spawn;


	//this->SpawnParticle();
	last_spawn = current_time;

	// total elapsed
	time_elapsed = current_time - last_loop;


	// -- new update method
	if (pHead->life > max_life)
	{
		memcpy(pHead, buffer, NUM_PARTICLES * sizeof(Particle));
	}
	Particle* t = pHead;
	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		t->Update(time_elapsed);
		t = t + 1;
	}

	last_loop = current_time;
}





void ParticleEmitter::draw()
{

	Particle* p = pHead;

	// move outside of the loop
	//Matrix tmp;
	Vect4D tmpPos;
	float cos;   
	float sin;    // 1.78


	for(int i = 0; i < NUM_PARTICLES; i++)
	{
		// --- new one
		//Particle local = *p;

		// get the position from this matrix



		// particle position
		tmpPos = p->position * 0.35f;
		cos = cosf(p->rotation);   // 2.74
		sin = sinf(p->rotation);    // 1.78

		/*Matrix tmp(_mm_set_ps(0.0f, 0.0f, (-sin) * p->scale.x * p->scale.y, cos * p->scale.x * p->scale.x),
			_mm_set_ps(0.0f, 0.0f, cos * p->scale.y * p->scale.x, sin * p->scale.y * p->scale.y),
			_mm_set_ps(0.0f, p->scale.z * p->scale.z, 0.0f, 0.0f),
			_mm_set_ps(1.0f, (inverseCameraMatrix.m14 + tmpPos.z) * p->scale.z, ((inverseCameraMatrix.m12 + tmpPos.x) * (-sin) + (inverseCameraMatrix.m13 + tmpPos.y) * cos) * p->scale.y, ((inverseCameraMatrix.m12 + tmpPos.x) * cos + (inverseCameraMatrix.m13 + tmpPos.y) * sin) * p->scale.x));*/
		
		Matrix tmp(cos * p->scale.x * p->scale.x, (-sin) * p->scale.x * p->scale.y, 0.0f, 0.0f,
			sin * p->scale.y * p->scale.y, cos * p->scale.y * p->scale.x, 0.0f, 0.0f,
			0.0f, 0.0f, p->scale.z * p->scale.z, 0.0f,
			((inverseCameraMatrix.m12 + tmpPos.x) * cos + (inverseCameraMatrix.m13 + tmpPos.y) * sin) * p->scale.x, ((inverseCameraMatrix.m12 + tmpPos.x) * (-sin) + (inverseCameraMatrix.m13 + tmpPos.y) * cos) * p->scale.y, (inverseCameraMatrix.m14 + tmpPos.z) * p->scale.z, 1.0f);
		
		glLoadMatrixf(reinterpret_cast<float*>(&(tmp)));

	//	// 1 -- tmp from BigFunction
	//								// vec		mat			  vec     float
		//matrixTransform(p->scale, cameraMatrix, tmpPos, p->rotation, tmp);
		//glLoadMatrixf(reinterpret_cast<float*>(&(tmp)));


		// draw the trangle strip
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);


		// increment to next point
		p += 1 ;
	}
}

void ParticleEmitter::Execute(Vect4D& pos, Vect4D& vel, Vect4D& sc)
{
	// Add some randomness...

	// --------------------------------------------------------------
	//   Ses it's ugly - I didn't write this so don't bitch at me   |
	//   Sometimes code like this is inside real commerical code    |
	//   ( so now you know how it feels )  |
	//---------------------------------------------------------------

	// x - variance
	float var = static_cast<float>(rand() % 1000) * 0.005f; // Var
	float sign = static_cast<float>(rand() % 2);  // Sign 
	float *t_pos = reinterpret_cast<float*>(&pos);
	float *t_var = &pos_variance[x];
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// y - variance
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// z - variance
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);

	// x  - add velocity
	t_pos = &vel[x];
	t_var = &vel_variance[x];
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// y - add velocity
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// z - add velocity
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	// correct the sign
	var = 1.5f * static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);

	if (sign == 0)
	{
		var *= -1.0;
	}
	sc = sc * var;
}


//Matrix ParticleEmitter::getFinalMatrix(const Matrix& scaleMatrix, const Matrix& transCamera, const Matrix& transParticle, const Matrix& rotParticle)
//{
//	// from Andre's slide WIMS Matrix multiplier
//	/*
//
//	// 		tmp = scaleMatrix * transCamera * transParticle *rotParticle * scaleMatrix;
//
//
//	scaleMatrix
//		{	sx		0		0		0 }
//		{	0		sy		0		0 }
//		{	0		0		sz		0 }
//		{	0		0		0		1 }
//
//	transCamera
//		{	1	0	0	0	}
//		{	0	1	0	0	}
//		{	0	0	1	0	}
//		{	cx	cy	cz	1	}
//
//		----------------------
//
//	scaleMatrix * transCamera
//		{   sx		0		0		0	}
//		{	0		sy		0		0	}
//		{	0		0		sz		0	}
//		{	cx		cy		cz		1	}
//
//	--------------------------------------------
//
//	transParticle
//		{	1	0	0	0	}
//		{	0	1	0	0	}
//		{	0	0	1	0	}
//		{	px	py	pz	1	}
//
//		----------------------
//
//	scaleMatrix * transCamera * transParticle
//
//		{   sx		0		0		0	}
//		{	0		sy		0		0	}
//		{	0		0		sz		0	}
//		{	cx+px	cy+py	cz+pz	1	}
//
//	--------------------------------------------
//
//	rotation matrix (rotParticle)
//		{	cos		-sin	0		0 }
//		{	sin		cos		0		0 }
//		{	0		0		1		0 }
//		{	0		0		0		1 }
//
//		----------------------
//
//	scaleMatrix * transCamera * transParticle * rotParticle
//
//		{   cos*sx							-sin*sx							0		0	}
//		{	ins*sy							cos*sy							0		0	}
//		{	0								0								sz		0	}
//		{	(cx+px)*cos+(cy+py)*sin		(cy+py)*(-sin)+(cy+py)*cos		cz+pz	1	}
//
//
//	--------------------------------------------
//	scaleMatrix
//		{	sx		0		0		0 }
//		{	0		sy		0		0 }
//		{	0		0		sz		0 }
//		{	0		0		0		1 }
//
//	tmp = scaleMatrix * transCamera * transParticle *rotParticle * scaleMatrix;
//
//		{   cos*sx*sx							-sin*sx*sy							0			0	}
//		{	sin*sy*sx							cos*sy*sy							0			0	}
//		{	0									0									sz*sz		0	}
//		{	((cx+px)*cos+(cy+py)*sin)*sx		((cx+px)*(-sin)+(cy+py)*cos)*sy		(cz+pz)*sz	1	}
//
//	sx = scaleMatrix.m0;
//	sy = scaleMatrix.m5;
//	sz = scaleMatrix.m10;
//
//	cos = rotParticle.m0;
//	-sin = rotParticle.m1;
//	sin = rotParticle.m4;
//
//	cx = transCamera.m12;
//	cy = transCamera.m13;
//	cz = transCamera.m14;
//
//	px = transParticle.m12;
//	py = transParticle.m13;
//	pz = transParticle.m14;
//
//	*/
//
//	Matrix tmp;
//	tmp.m0 = scaleMatrix.m0 * scaleMatrix.m0 * rotParticle.m0;
//	tmp.m1 = rotParticle.m1 * scaleMatrix.m0 * scaleMatrix.m5;
//	tmp.m2 = 0.0f;
//	tmp.m3 = 0.0f;
//	tmp.m4 = rotParticle.m4 * scaleMatrix.m5 * scaleMatrix.m0;
//	tmp.m5 = rotParticle.m0 * scaleMatrix.m5 * scaleMatrix.m5;
//	tmp.m6 = 0.0f;
//	tmp.m7 = 0.0f;
//	tmp.m8 = 0.0f;
//	tmp.m9 = 0.0f;
//	tmp.m10 = scaleMatrix.m10 * scaleMatrix.m10;
//	tmp.m11 = 0.0f;
//	tmp.m12 = ((transCamera.m12 + transParticle.m12) * rotParticle.m0
//		+ (transCamera.m13 + transParticle.m13) * rotParticle.m4)
//		* scaleMatrix.m0;
//	tmp.m13 = ((transCamera.m12 + transParticle.m12) * rotParticle.m1
//		+ (transCamera.m13 + transParticle.m13) * rotParticle.m0)
//		* scaleMatrix.m5;
//	tmp.m14 = (transCamera.m14 + transParticle.m14) * scaleMatrix.m10;
//	tmp.m15 = 1.0f;
//	return tmp;
//}
//
//void ParticleEmitter::matrixTransform(const Vect4D& scale, const Matrix& transCam, const Vect4D& transPar, const float& rot, Matrix& out)
//{
//	float cos = cosf(rot);   // 2.74
//	float minusSin = -sinf(rot);    // 2.07
//	float sin = sinf(rot);    // 1.78
//	/*
//	tmp = scaleMatrix * transCamera * transParticle * rotParticle * scaleMatrix;
//
//	{   cos* sx* sx									-sin * sx * sy									0				0	}
//	{	sin* sy* sx									cos* sy* sy										0				0	}
//	{	0											0												sz * sz			0	}
//	{	((cx + px) * cos + (cy + py) * sin)* sx		((cx + px) * (-sin) + (cy + py) * cos)* sy		(cz + pz)* sz	1	}
//	*/
//
//	out.m0 = cos * scale.x * scale.x;
//	out.m1 = minusSin * scale.x * scale.y;
//	out.m2 = 0.0f;
//	out.m3 = 0.0f;
//	out.m4 = sin * scale.y * scale.x;
//	out.m5 = cos * scale.y * scale.y;
//	out.m6 = 0.0f;
//	out.m7 = 0.0f;
//	out.m8 = 0.0f;
//	out.m9 = 0.0f;
//	out.m10 = scale.z * scale.z;
//	out.m11 = 0.0f;
//	out.m12 = ((transCam.m12 + transPar.x) * cos + (transCam.m13 + transPar.y) * sin) * scale.x;
//	out.m13 = ((transCam.m12 + transPar.x) * minusSin + (transCam.m13 + transPar.y) * cos) * scale.y;
//	out.m14 = (transCam.m14 + transPar.z) * scale.z;
//	out.m15 = 1.0f;
//
//	// ref for the big transform
//	// 2 -- pass everything through matrix c-tor
//	//	/*tmp.m0 = cos * p->scale.x * p->scale.x;
//	//	tmp.m1 = minusSin * p->scale.x * p->scale.y;
//	//	tmp.m2 = 0.0f;
//	//	tmp.m3 = 0.0f;
//
//	//	tmp.m4 = sin * p->scale.y * p->scale.x;
//	//	tmp.m5 = cos * p->scale.y * p->scale.y;
//	//	tmp.m6 = 0.0f;
//	//	tmp.m7 = 0.0f;
//
//	//	tmp.m8 = 0.0f;
//	//	tmp.m9 = 0.0f;
//	//	tmp.m10 = p->scale.z * p->scale.z;
//	//	tmp.m11 = 0.0f;
//
//	//	tmp.m12 = ((cameraMatrix.m12 + tmpPos.x) * cos + (cameraMatrix.m13 + tmpPos.y) * sin) * p->scale.x;
//	//	tmp.m13 = ((cameraMatrix.m12 + tmpPos.x) * minusSin + (cameraMatrix.m13 + tmpPos.y) * cos) * p->scale.y;
//	//	tmp.m14 = (cameraMatrix.m14 + tmpPos.z) * p->scale.z;
//	//	tmp.m15 = 1.0f;*/
//}

// End of file
