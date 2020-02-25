//----------------------------------------------------------------------------
// Copyright 2019, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

// This is a 4 dimensional Vect4D class
#ifndef Vect4D_H
#define Vect4D_H

// includes
#include "Enum.h"
#include <xmmintrin.h>
#include <smmintrin.h>  

// Foward Declarations
class Matrix;

// class
class Vect4D: Align16
{
public:
	friend class Matrix;

	// big four
	Vect4D();
	//Vect4D(const Vect4D&) = default;
	~Vect4D() = default;
	//Vect4D& operator= (const Vect4D& tmp) = default;
	//Vect4D& operator= (const Vect4D& t);
	Vect4D& operator= (const __m128& tmp);

	Vect4D(const __m128& tm) : _m(tm) {};

	//Vect4D(float tx, float ty, float tz, float tw = 1.0f) : x(tx), y(ty), z(tz), w(tw) {};
	Vect4D(float tx, float ty, float tz, float tw = 1.0f);



	void set(const float& tx, const float& ty, const float& tz, const float tw = 1.0f);

	const Vect4D operator + (const Vect4D& t) const;
	const Vect4D operator - (const Vect4D& t)const;
	const Vect4D operator * (const float& scale)  ;
	//Vect4D& operator += (Vect4D rhs);
	//const Vect4D& Vect4D::operator*=(const float& scale);

	//const Vect4D norm(Vect4D& out) const;

	const void norm(Vect4D& out)const ;
	const void Cross(const Vect4D &vin, Vect4D &vout) const;

	float &operator[](const VECT_ENUM e) ;

public:

	// anonymous union
	union
	{
		__m128	_m;

		// anonymous struct
		struct
		{
			float x;
			float y;
			float z;
			float w;
		};
	};
};


#endif  

// End of file

