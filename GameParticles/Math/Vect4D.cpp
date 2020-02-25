//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include <math.h>
#include "Vect4D.h"

Vect4D::Vect4D()
{
	this->x = 0.0f;
	this->y = 0.0f;
	this->z = 0.0f;
	this->w = 1.0f;
	//this->_m = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);

}

Vect4D::Vect4D(float tx, float ty, float tz, float tw)
	//: x(tx), y(ty), z(tz), w(tw)
{
	// amo commented out these
	this->x = tx;
	this->y = ty;
	this->z = tz;
	this->w = tw;
	
}



//// assginment operator overloading: input: Vect4D
//Vect4D& Vect4D::operator= (const Vect4D& tmp) 
//{
//	this->_m = _mm_set_ps(tmp.w, tmp.z, tmp.y, tmp.x); //4.99
//	return *this;
//}

Vect4D& Vect4D::operator = (const __m128& tmp)
{

	// ver 1
	//this->_m = tmp;

	// ver 2
	//this->w = tmp.m128_f32[0];
	//this->z = tmp.m128_f32[1];
	//this->y = tmp.m128_f32[2];
	//this->x = tmp.m128_f32[3];

	// ver 3
	this->_m.m128_f32[0] = tmp.m128_f32[0];
	this->_m.m128_f32[1] = tmp.m128_f32[1];
	this->_m.m128_f32[2] = tmp.m128_f32[2];
	this->_m.m128_f32[3] = tmp.m128_f32[3];


	return *this;
}

//Vect4D& Vect4D::operator = (const Vect4D& t)
//{
//	this->x = t.x;
//	this->y = t.y;
//	this->z = t.z;
//	this->w = t.w;
//	return *this;
//}




const void Vect4D::norm(Vect4D& out) const
{
	// amo:: changed from sqrt to sqrtf
	float mag = sqrtf(this->x * this->x + this->y * this->y + this->z * this->z);
	mag = 1.0f / mag;

	// original
	//out.x = this->x * mag;
	//out.y = this->y * mag;
	//out.z = this->z * mag;
	//out.w = 1.0f;


	// SIMD ver.
	out._m = _mm_setr_ps(this->x * mag, this->y * mag, this->z * mag, 1.0f);

	// amo: rm this part, wouldn't happen
	//if (0.0f < mag)
	//{
	//	
	//}
}

const Vect4D Vect4D::operator + (const Vect4D& t) const
{
	// 0. original
	//Vect4D out;

	//out.x = this->x + t.x;
	//out.y = this->y + t.y;
	//out.z = this->z + t.z;

	//return out;

	// 1. RVO
	return Vect4D(this->x + t.x, this->y + t.y, this->z + t.z, 1.0f);

	// 2. SIMD
	//__m128 mt = _mm_set_ps(t.w, t.z, t.y, t.x);
	//return Vect4D(_mm_add_ps(this->_m, mt));

}

const Vect4D Vect4D::operator - (const Vect4D& t) const
{
	//Vect4D out;

	//out.x = this->x - t.x;
	//out.y = this->y - t.y;
	//out.z = this->z - t.z;

	//return out;


	//Vect4D tmp;

	//__m128 mt = _mm_set_ps(t.w, t.z, t.y, t.x);
	//__m128 minus = _mm_sub_ps(this->_m, mt);
	//tmp._m = minus;

	//return tmp;

	// 2 --- RVO version (faster)
	//return Vect4D(this->x - t.x, this->y - t.y, this->z - t.z, 1.0f);

	// 3 --- SIMD + RVO
	__m128 mt = _mm_set_ps(t.w, t.z, t.y, t.x);
	return Vect4D(_mm_sub_ps(this->_m, mt));

}

const Vect4D Vect4D::operator *(const float& scale) 
{
	// 1. --- original ver.
	//Vect4D tmp;

	//tmp.x = this->x * scale;
	//tmp.y = this->y * scale;
	//tmp.z = this->z * scale;

	//return tmp;

	// 2. --- RVO ver.
	return Vect4D(this->x * scale, this->y * scale, this->z * scale, 1.0f);  //140 5. 52% --> 3.89%

	// 3 --- call =operator(__m128) ver.
	//Vect4D tmp;   // 1.56% -> 1.39%

	//__m128 _scale = _mm_set_ps1(scale);
	//tmp._m = _mm_mul_ps(this->_m, _scale);


	//__m128 _scale = _mm_set_ps1(scale);
	//__m128 mul = _mm_mul_ps(this->_m, _scale);
	//this->_m = mul;
	//tmp._m = mul;

	// 4 ---
	//__m128 _scale = _mm_set_ps(1.0f,scale, scale, scale);
	//return Vect4D(_mm_mul_ps(this->_m, _scale));   	// total: 4%

	//this->_m = _mm_mul_ps(this->_m, _scale);
	//return *this;

}

//Vect4D& Vect4D::operator+=(Vect4D rhs)
//{
//	this->_m = (_mm_add_ps(this->_m, rhs._m));
//	return *this;
//}
//
//const Vect4D& Vect4D::operator*=(const float& scale)
//{
//	*this = Vect4D(this->x * scale, this->y * scale, this->z * scale, 0.0f);
//	return *this;
//}

float& Vect4D::operator[](const VECT_ENUM e) 
{
	switch (e)
	{
	case 0:
		return x;
		break;
	case 1:
		return y;
		break;
	case 2:
		return z;
		break;
	case 3:
		return w;
		break;
	default:
		assert(0);
		return x;
	}
}

const void Vect4D::Cross(const Vect4D& vin, Vect4D& vout) const
{
	// original ver.
	vout.x = (y*vin.z - z * vin.y);
	vout.y = (z*vin.x - x * vin.z);
	vout.z = (x*vin.y - y * vin.x);
	//vout.w = 1.0f;

	// SIMD ver.

	//__m128 pos1 = _mm_set_ps(1.0f, this->x, this->z, this->y);
	//__m128 vin1 = _mm_set_ps(1.0f, vin.y, vin.x, vin.z);
	//__m128 pos2 = _mm_set_ps(1.0f, this->y, this->x, this->z);
	//__m128 vin2 = _mm_set_ps(1.0f, vin.x, vin.z, vin.y);

	//__m128 mul1 = _mm_mul_ps(pos1, vin1);  // 0.24
	//__m128 mul2 = _mm_mul_ps(pos2, vin2); // 0.45
	//vout._m = _mm_sub_ps(mul1, mul2); // 1.09
	//vout.w = 1.0f;  // 0,1
}

void Vect4D::set(const float& tx, const float& ty, const float& tz, const float tw)
{
	// original
	//this->x = tx;
	//this->y = ty;
	//this->z = tz;
	//this->w = tw;

	// SIMD
	this->_m = _mm_set_ps(tw, tz, ty, tx);
}

// End of file
