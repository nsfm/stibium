#ifndef MATH_F_H
#define MATH_F_H

#include <math.h>

#include "fab/tree/math/math_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @file tree/math/math_f.h
    @brief Functions for doing math on floating-point numbers.
    @details These functions take in input floats A and B,
    and return the result of their computation.
*/

// Binary functions
inline float add_f(float A, float B) { return A+B; }
inline float sub_f(float A, float B) { return A-B; }
inline float mul_f(float A, float B) { return A*B; }
inline float div_f(float A, float B) { return A/B; }
/*  Bit-exact replicas of this platform's libm fmin/fmax (x86 minss
 *  semantics: ties including ±0 return B, NaN-vs-NaN returns A),
 *  minus sNaN quieting, which arithmetic can never produce.  Written
 *  out so they inline and vectorize instead of a per-voxel PLT call
 *  into libm — min/max chains are most of a CSG model's tape.  */
inline float min_f(float A, float B)
{
    if (isnan(B)) return A;
    if (isnan(A)) return B;
    return A < B ? A : B;
}
inline float max_f(float A, float B)
{
    if (isnan(B)) return A;
    if (isnan(A)) return B;
    return A > B ? A : B;
}
inline float pow_f(float A, float B) { return pow(A, B); }

////////////////////////////////////////////////////////////////////////////////

inline float abs_f(float A) { return fabs(A); }
inline float square_f(float A) { return A*A; }
inline float sqrt_f(float A) { return A < 0 ? 0 : sqrt(A); }
inline float sin_f(float A) { return sin(A); }
inline float cos_f(float A) { return cos(A); }
inline float tan_f(float A) { return tan(A); }

inline float asin_f(float A)
{
    if (A < -1)     return -M_PI_2;
    else if (A > 1) return M_PI_2;
    else            return asin(A);
}

inline float acos_f(float A)
{
    if (A < -1)     return M_PI;
    else if (A > 1) return 0;
    else            return acos(A);
}

inline float atan_f(float A) { return atan(A); }
inline float atan2_f(float A, float B) { return atan2(A, B); }

// GLSL-style mod: result has the sign of B, mod(A, B) in [0, B) for B > 0
inline float mod_f(float A, float B)
{
    if (B == 0)     return 0;
    return A - B * floor(A / B);
}
inline float floor_f(float A) { return floor(A); }

// Natural log, clamped at the smallest normal float so non-positive
// inputs yield a large negative number instead of -inf/NaN
inline float log_f(float A)
{
    return log(A > 1.17549e-38f ? A : 1.17549e-38f);
}
inline float neg_f(float A) { return -A; }
inline float exp_f(float A) { return exp(A); }

////////////////////////////////////////////////////////////////////////////////

inline float X_f(float X) { return X; }
inline float Y_f(float Y) { return Y; }
inline float Z_f(float Z) { return Z; }

#ifdef __cplusplus
}
#endif

#endif
