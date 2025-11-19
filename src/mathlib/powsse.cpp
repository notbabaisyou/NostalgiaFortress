//========= Copyright Valve Corporation, All rights reserved. ============//
//
// Purpose: 
//
//=====================================================================================//

#include "mathlib/ssemath.h"

// NOTE: This has to be the last file included!
#include "tier0/memdbgon.h"


fltx4 Pow_FixedPoint_Exponent_SIMD( const fltx4 & x, int exponent)
{
	fltx4 rslt=Four_Ones;									// x^0=1.0
	int xp=abs(exponent);
	if (xp & 3)												// fraction present?
	{
		fltx4 sq_rt=SqrtEstSIMD(x);
		if (xp & 1)											// .25?
			rslt=SqrtEstSIMD(sq_rt);						// x^.25
		if (xp & 2)
			rslt=MulSIMD(rslt,sq_rt);
	}
	xp>>=2;													// strip fraction
	fltx4 curpower=x;										// curpower iterates through  x,x^2,x^4,x^8,x^16...

	while(1)
	{
		if (xp & 1)
			rslt=MulSIMD(rslt,curpower);
		xp>>=1;
		if (xp)
			curpower=MulSIMD(curpower,curpower);
		else
			break;
	}
	if (exponent<0)
		return ReciprocalEstSaturateSIMD(rslt);				// pow(x,-b)=1/pow(x,b)
	else
		return rslt;
}

float FastLog2(float val)
{
	union { float val; int32_t x; } u = { val, };
    float log_2 = (float)(((u.x >> 23) & 255) - 128);              
    u.x   &= ~(255 << 23);
    u.x   += 127 << 23;
    log_2 += ((-0.34484843f) * u.val + 2.02466578f) * u.val - 0.67487759f; 
    return (log_2);
}

float FastPow2(float i)
{
	__m128 v = _mm_set_ss(i);
	v = _mm_mul_ss(v, v);
	return _mm_cvtss_f32(v);
}

float FastPow(float a, float b)
{
	return FastPow2(b * FastLog2(a));
}

float FastPow10( float i )
{
	return FastPow2( i * 3.321928f );
}

