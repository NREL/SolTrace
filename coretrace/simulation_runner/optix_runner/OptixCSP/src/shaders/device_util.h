#pragma once

#include <vector_types.h>
#include <vector_functions.h>  

#include <cmath>
#include <cstdlib>

#define INLINE __forceinline__
#define HOSTDEVICE __host__ __device__

#ifndef M_PIf
#define M_PIf       3.14159265358979323846f
#endif
#ifndef M_PI_2f
#define M_PI_2f     1.57079632679489661923f
#endif

INLINE HOSTDEVICE float3 make_float3(const float4& v0) { return make_float3(v0.x, v0.y, v0.z); }

/** float2 math **/
/** dot product **/
INLINE HOSTDEVICE float dot(const float2& a, const float2& b)
{
    return a.x * b.x + a.y * b.y;
}

/** float3 math **/
/** negate */
INLINE HOSTDEVICE float3 operator-(const float3& a)
{
    return make_float3(-a.x, -a.y, -a.z);
}

/** @} */
/** add
* @{
*/
INLINE HOSTDEVICE float3 operator+(const float3& a, const float3& b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
INLINE HOSTDEVICE float3 operator+(const float3& a, const float b)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}
INLINE HOSTDEVICE float3 operator+(const float a, const float3& b)
{
    return make_float3(a + b.x, a + b.y, a + b.z);
}
INLINE HOSTDEVICE void operator+=(float3& a, const float3& b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}
/** @} */

/** subtract
* @{
*/
INLINE HOSTDEVICE float3 operator-(const float3& a, const float3& b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
INLINE HOSTDEVICE float3 operator-(const float3& a, const float b)
{
    return make_float3(a.x - b, a.y - b, a.z - b);
}
INLINE HOSTDEVICE float3 operator-(const float a, const float3& b)
{
    return make_float3(a - b.x, a - b.y, a - b.z);
}
INLINE HOSTDEVICE void operator-=(float3& a, const float3& b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

/** multiply
* @{
*/
INLINE HOSTDEVICE float3 operator*(const float3& a, const float3& b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
INLINE HOSTDEVICE float3 operator*(const float3& a, const float s)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
INLINE HOSTDEVICE float3 operator*(const float s, const float3& a)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
INLINE HOSTDEVICE void operator*=(float3& a, const float3& s)
{
    a.x *= s.x; a.y *= s.y; a.z *= s.z;
}
INLINE HOSTDEVICE void operator*=(float3& a, const float s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

INLINE HOSTDEVICE float dot(const float3& a, const float3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

INLINE HOSTDEVICE float3 cross(const float3& a, const float3& b)
{
    return make_float3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

INLINE HOSTDEVICE float length(const float3& v)
{
    return sqrtf(dot(v, v));
}

INLINE HOSTDEVICE float3 normalize(const float3& v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

INLINE HOSTDEVICE float4 make_float4(const float3& v, float w)
{
    return make_float4(v.x, v.y, v.z, w);
}

INLINE HOSTDEVICE float4 make_float4(const float v0, const float3& v1) { return make_float4(v0, v1.x, v1.y, v1.z); }


/** Faceforward
* Returns N if dot(i, nref) > 0; else -N;
* Typical usage is N = faceforward(N, -ray.dir, N);
* Note that this is opposite of what faceforward does in Cg and GLSL */
INLINE HOSTDEVICE float3 faceforward(const float3& n, const float3& i, const float3& nref)
{
  // return n * copysignf(1.0f, dot(i, nref));
  return dot(i, nref) >= 0.0f ? n : -n;
}


#define float3_as_args( u )                                                                                            \
    reinterpret_cast<unsigned int&>( ( u ).x ), reinterpret_cast<unsigned int&>( ( u ).y ),                            \
        reinterpret_cast<unsigned int&>( ( u ).z )
