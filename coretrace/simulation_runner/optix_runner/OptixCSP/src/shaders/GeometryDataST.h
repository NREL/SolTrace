#pragma once
#include "device_util.h"
#include "soltrace_constants.h"
#ifndef __CUDACC_RTC__
#include <cassert>
#include <cstdint>
#else
#define assert(x) /*nop*/
#endif

// TODO: get rid of ST suffix, no clue what it was for ...

namespace OptixCSP
{
    struct GeometryDataST
    {
        enum Type
        {
            PARALLELOGRAM = 0,
            CYLINDER_Y = 1,
            RECTANGLE_PARABOLIC = 2,
            UNKNOWN_TYPE = 3,
            RECTANGLE_FLAT = 4,
            TRIANGLE_FLAT = 5,
            QUADRILATERAL_FLAT = 6,
            CIRCLE_FLAT = 7,
            HEXAGON_FLAT = 8,
            ANNULUS_FLAT = 9,
            CIRCLE_PARABOLIC = 10,
            HEXAGON_PARABOLIC = 11,
            TRIANGLE_PARABOLIC = 12,
            ANNULUS_PARABOLIC = 13,
            QUADRILATERAL_PARABOLIC = 14,
            RECTANGLE_SPHERICAL = 15,
            CIRCLE_SPHERICAL = 16,
            HEXAGON_SPHERICAL = 17,
            ANNULUS_SPHERICAL = 18,
            TRIANGLE_SPHERICAL = 19,
            QUADRILATERAL_SPHERICAL = 20
        };

        struct Parallelogram
        {
            Parallelogram() = default;
            Parallelogram(float3 v1, float3 v2, float3 anchor)
                : v1(v1), v2(v2), anchor(anchor)
            {
                float3 normal = normalize(cross(v1, v2));
                float d = dot(normal, anchor);
                this->v1 *= 1.0f / dot(v1, v1);
                this->v2 *= 1.0f / dot(v2, v2);
                plane = make_float4(normal, d);
            }
            float4 plane;  // plane equation: (normal, dot(anchor, normal))
            float3 v1;     // edge vector 1, stored as v1/dot(v1,v1)
            float3 v2;     // edge vector 2, stored as v2/dot(v2,v2)
            float3 anchor; // corner point of the parallelogram
        };

        // same as parallelogram, however defined with different attributes
        struct Rectangle_Flat
        {
            Rectangle_Flat() = default;
            Rectangle_Flat(float3 center, float3 x, float3 y, float width, float height)
                : center(center), x(x), y(y), width(width), height(height)
            {
                float3 normal = normalize(cross(x, y));
                float d = dot(normal, center);
                plane = make_float4(normal, d);
            }

            float4 plane;  // normal unit vector, dot(center, normal)
            float3 center; // center in global coordinates
            float3 x;      // local x axis unit vector
            float3 y;      // local y axis unit vector
            float width;   // full width along x
            float height;  // full height along y
        };

        struct Cylinder_Y
        {
            Cylinder_Y() = default;
            Cylinder_Y(float3 center, float radius, float half_height, float3 base_x, float3 base_z)
                : center(center), radius(radius), half_height(half_height), base_x(base_x), base_z(base_z)
            {
                assert(dot(base_x, base_z) < 1e-3f);
            }

            float3 center;     // center of the cylinder in global coordinates
            float radius;      // radius
            float half_height; // half the height along the local y axis
            float3 base_x;     // x axis of the cylinder
            float3 base_z;     // z axis of the cylinder
        };

        struct Rectangle_Parabolic
        {

            Rectangle_Parabolic() = default;
            Rectangle_Parabolic(float3 v1, float3 v2, float3 anchor, float curv_x, float curv_y)
                : v1(v1), v2(v2), anchor(anchor), curv_x(curv_x), curv_y(curv_y)
            {
                float3 normal = normalize(cross(v1, v2));
                float d = dot(normal, anchor);
                this->v1 *= 1.0f / dot(v1, v1);
                this->v2 *= 1.0f / dot(v2, v2);
                plane = make_float4(normal, d);
            }

            float4 plane;  // plane equation of the base rectangle: (normal, dot(anchor, normal))
            float3 v1;     // edge vector 1, stored as v1/dot(v1,v1)
            float3 v2;     // edge vector 2, stored as v2/dot(v2,v2)
            float3 anchor; // corner point of the base rectangle
            // float3 focus;
            float curv_x; // curvature along local x axis
            float curv_y; // curvature along local y axis
        };

        struct Triangle_Flat
        {
            Triangle_Flat() = default;
            Triangle_Flat(const float3 &a, const float3 &b, const float3 &c)
                : v0(a), e1(b - a), e2(c - a)
            {
                normal = normalize(cross(e1, e2));
                d = dot(normal, v0);
            }
            float3 v0;     // base vertex
            float3 e1, e2; // edges
            float3 normal; // normal unit vector
            float d;       // plane distance
        };

        struct Quadrilateral_Flat
        {
            Quadrilateral_Flat() = default;
            Quadrilateral_Flat(const float3 &a, const float3 &b,
                               const float3 &c, const float3 &d)
                : p0(a), p1(b), p2(c), p3(d)
            {
                float3 e1 = p1 - p0;
                float3 e2 = p3 - p0;
                normal = normalize(cross(e1, e2));
            }
            float3 p0, p1, p2, p3; // Vertices in counterclockwise order
            float3 normal;         // Positive direction follows right-hand rule
        };

        struct Circle_Flat
        {
            Circle_Flat() = default;
            // Circle_Flat(const float radius) : r(radius) {}
            Circle_Flat(const float3 &origin, const float3 &normal, const float &radius)
                : r(radius), center(origin)
            {
                const float3 n = normalize(normal);
                plane = make_float4(n, dot(center, n));
            }
            float4 plane;  // normal unit vector, dot(center, normal)
            float3 center; // local origin in global coordinates
            float r;       // radius
        };

        struct Hexagon_Flat
        {
            Hexagon_Flat() = default;
            Hexagon_Flat(const float3 &origin, const float3 &normal,
                         const float3 &x_ax, const float3 &y_ax,
                         const float &side_length)
                : center(origin), x_axis(x_ax), y_axis(y_ax), s(side_length)
            {
                const float3 n = normalize(normal);
                plane = make_float4(n, dot(center, n));
            }
            float4 plane;  // normal unit vector, dot(center, normal)
            float3 center; // local origin in global coordinates
            float3 x_axis; // unit vector
            float3 y_axis; // unit vector
            float s;       // side length
        };

        struct Annulus_Flat
        {
            Annulus_Flat() = default;
            Annulus_Flat(const float3 &origin, const float3 &normal,
                         const float3 &x_ax, const float3 &y_ax,
                         const float &r_inner, const float &r_outer,
                         const float &arc)
                : center(origin), x_axis(x_ax), y_axis(y_ax),
                  ri(r_inner), ro(r_outer), arc(arc)
            {
                const float3 n = normalize(normal);
                plane = make_float4(n, dot(center, n));
            }
            float4 plane;  // normal unit vector, dot(center, normal)
            float3 center; // local origin in global coordinates
            float3 x_axis; // local x axis unit vector (arc is centered about this axis)
            float3 y_axis; // local y axis unit vector
            float ri;      // inner radius
            float ro;      // outer radius
            float arc;     // total arc angle in radians, centered on x_axis
        };

        struct Circle_Parabolic
        {
            Circle_Parabolic() = default;
            Circle_Parabolic(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                             const float &curv_x, const float &curv_y, const float &r)
                : center(origin), x_axis(x_ax), y_axis(y_ax),
                  cx(curv_x), cy(curv_y), radius(r)
            {
            }
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float cx;
            float cy;
            float radius;
        };

        struct Hexagon_Parabolic
        {
            Hexagon_Parabolic() = default;
            Hexagon_Parabolic(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                              const float &curv_x, const float &curv_y, const float &side_len)
                : center(origin), x_axis(x_ax), y_axis(y_ax),
                  cx(curv_x), cy(curv_y), s(side_len)
            {
            }
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float cx;
            float cy;
            float s;
        };

        struct Triangle_Parabolic
        {
            Triangle_Parabolic() = default;
            // Vertices v0, v1, v2 are in the local XY aperture frame.
            // The constructor precomputes the barycentric inverse transform so
            // the aperture test reduces to two dot products with no per-ray division.
            Triangle_Parabolic(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                               const float &curv_x, const float &curv_y,
                               const float2 &v0, const float2 &v1, const float2 &v2)
                : center(origin), x_axis(x_ax), y_axis(y_ax),
                  cx(curv_x), cy(curv_y)
            {
                const float2 e1 = make_float2(v1.x - v0.x, v1.y - v0.y);
                const float2 e2 = make_float2(v2.x - v0.x, v2.y - v0.y);
                const float inv_det = 1.0f / (e1.x * e2.y - e1.y * e2.x);
                // u = dot(utest, float3(px, py, 1.0f))
                utest = make_float3(e2.y, -e2.x, v0.y * e2.x - v0.x * e2.y) * inv_det;
                // v = dot(vtest, float3(px, py, 1.0f))
                vtest = make_float3(-e1.y, e1.x, v0.x * e1.y - v0.y * e1.x) * inv_det;
            }
            float3 center; // element origin in global coordinates
            float3 x_axis; // local x axis unit vector
            float3 y_axis; // local y axis unit vector
            float cx;      // curvature along local x axis
            float cy;      // curvature along local y axis
            float3 utest;  // precomputed row: u = dot(utest, float3(px, py, 1.0f))
            float3 vtest;  // precomputed row: v = dot(vtest, float3(px, py, 1.0f))
        };

        struct Annulus_Parabolic
        {
            Annulus_Parabolic() = default;
            Annulus_Parabolic(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                              const float &curv_x, const float &curv_y,
                              const float &r_inner, const float &r_outer, const float &arc)
                : center(origin), x_axis(x_ax), y_axis(y_ax),
                  cx(curv_x), cy(curv_y),
                  ri(r_inner), ro(r_outer), arc(arc)
            {
            }
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float cx;
            float cy;
            float ri;
            float ro;
            float arc; // in radians
        };

        struct Quadrilateral_Parabolic
        {
            Quadrilateral_Parabolic() = default;
            Quadrilateral_Parabolic(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                                    const float &curv_x, const float &curv_y,
                                    const float2 &v0, const float2 &v1,
                                    const float2 &v2, const float2 &v3)
                : center(origin), x_axis(x_ax), y_axis(y_ax), cx(curv_x), cy(curv_y)
            {
                const float2 e1 = make_float2(v1.x - v0.x, v1.y - v0.y);
                const float2 e2 = make_float2(v2.x - v0.x, v2.y - v0.y);
                const float inv_det = 1.0f / (e1.x * e2.y - e1.y * e2.x);
                // u = dot(utest, float3(px, py, 1.0f))
                u1test = make_float3(e2.y, -e2.x, v0.y * e2.x - v0.x * e2.y) * inv_det;
                // v = dot(vtest, float3(px, py, 1.0f))
                v1test = make_float3(-e1.y, e1.x, v0.x * e1.y - v0.y * e1.x) * inv_det;

                const float2 e3 = make_float2(v3.x - v0.x, v3.y - v0.y);
                const float inv_det2 = 1.0f / (e2.x * e3.y - e2.y * e3.x);
                u2test = make_float3(e3.y, -e3.x, v0.y * e3.x - v0.x * e3.y) * inv_det2;
                v2test = make_float3(-e2.y, e2.x, v0.x * e2.y - v0.y * e2.x) * inv_det2;
            }
            float3 center; // element origin in global coordinates
            float3 x_axis; // local x axis unit vector
            float3 y_axis; // local y axis unit vector
            float cx;      // curvature along local x axis
            float cy;      // curvature along local y axis
            float3 u1test;  // precomputed row: u = dot(utest, float3(px, py, 1.0f))
            float3 v1test;  // precomputed row: v = dot(vtest, float3(px, py, 1.0f))
            float3 u2test;
            float3 v2test;
        };

        // -----------------------------------------------------------------------
        // Spherical surface data structures.
        // All spherical aperture structs share the sphere radius R.
        // The surface equation in the local element frame is:
        //   z(x, y) = (x^2 + y^2) / [(R + sqrt(R^2 - (x^2 + y^2)))]
        // which is the lower cap of a sphere with radius R centred at
        // (0, 0, R) in element-local coordinates.
        // -----------------------------------------------------------------------

        struct Rectangle_Spherical
        {
            Rectangle_Spherical() = default;
            Rectangle_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                                 const float &radius, const float &w, const float &h,
                                 const float &xc, const float &yc)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(radius),
                  width(w), height(h), x_coord(xc), y_coord(yc)
            {}
            float3 center;  // element origin in global coordinates
            float3 x_axis;  // local x axis unit vector
            float3 y_axis;  // local y axis unit vector
            float R;        // sphere radius
            float width;    // full width along x
            float height;   // full height along y
            float x_coord;  // aperture x offset
            float y_coord;  // aperture y offset
        };

        struct Circle_Spherical
        {
            Circle_Spherical() = default;
            Circle_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                             const float &sphere_R, const float &r)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(sphere_R), radius(r)
            {}
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float R;
            float radius;
        };

        struct Hexagon_Spherical
        {
            Hexagon_Spherical() = default;
            Hexagon_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                              const float &radius, const float &side_len)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(radius), s(side_len)
            {}
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float R;
            float s;  // circumradius (vertex-to-center distance)
        };

        struct Annulus_Spherical
        {
            Annulus_Spherical() = default;
            Annulus_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                              const float &radius,
                              const float &r_inner, const float &r_outer, const float &arc_ang)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(radius),
                  ri(r_inner), ro(r_outer), arc(arc_ang)
            {}
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float R;
            float ri;
            float ro;
            float arc;  // in radians
        };

        struct Triangle_Spherical
        {
            Triangle_Spherical() = default;
            Triangle_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                               const float &radius,
                               const float2 &v0, const float2 &v1, const float2 &v2)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(radius)
            {
                const float2 e1 = make_float2(v1.x - v0.x, v1.y - v0.y);
                const float2 e2 = make_float2(v2.x - v0.x, v2.y - v0.y);
                const float inv_det = 1.0f / (e1.x * e2.y - e1.y * e2.x);
                utest = make_float3(e2.y, -e2.x, v0.y * e2.x - v0.x * e2.y) * inv_det;
                vtest = make_float3(-e1.y, e1.x, v0.x * e1.y - v0.y * e1.x) * inv_det;
            }
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float R;
            float3 utest;
            float3 vtest;
        };

        struct Quadrilateral_Spherical
        {
            Quadrilateral_Spherical() = default;
            Quadrilateral_Spherical(const float3 &origin, const float3 &x_ax, const float3 &y_ax,
                                    const float &radius,
                                    const float2 &v0, const float2 &v1,
                                    const float2 &v2, const float2 &v3)
                : center(origin), x_axis(x_ax), y_axis(y_ax), R(radius)
            {
                const float2 e1 = make_float2(v1.x - v0.x, v1.y - v0.y);
                const float2 e2 = make_float2(v2.x - v0.x, v2.y - v0.y);
                const float inv_det = 1.0f / (e1.x * e2.y - e1.y * e2.x);
                u1test = make_float3(e2.y, -e2.x, v0.y * e2.x - v0.x * e2.y) * inv_det;
                v1test = make_float3(-e1.y, e1.x, v0.x * e1.y - v0.y * e1.x) * inv_det;

                const float2 e3 = make_float2(v3.x - v0.x, v3.y - v0.y);
                const float inv_det2 = 1.0f / (e2.x * e3.y - e2.y * e3.x);
                u2test = make_float3(e3.y, -e3.x, v0.y * e3.x - v0.x * e3.y) * inv_det2;
                v2test = make_float3(-e2.y, e2.x, v0.x * e2.y - v0.y * e2.x) * inv_det2;
            }
            float3 center;
            float3 x_axis;
            float3 y_axis;
            float R;
            float3 u1test;
            float3 v1test;
            float3 u2test;
            float3 v2test;
        };

        GeometryDataST() = default;

        void setParallelogram(const Parallelogram &p)
        {
            assert(type == UNKNOWN_TYPE);
            type = PARALLELOGRAM;
            parallelogram = p;
        }

        __host__ __device__ const Parallelogram &getParallelogram() const
        {
            assert(type == PARALLELOGRAM);
            return parallelogram;
        }

        void setRectangle_Flat(const Rectangle_Flat &r)
        {
            assert(type == UNKNOWN_TYPE);
            type = RECTANGLE_FLAT;
            rectangle_flat = r;
        }

        __host__ __device__ const Rectangle_Flat &getRectangle_Flat() const
        {
            assert(type == RECTANGLE_FLAT);
            return rectangle_flat;
        }

        void setCylinder_Y(const Cylinder_Y &c)
        {
            assert(type == UNKNOWN_TYPE);
            type = CYLINDER_Y;
            cylinder_y = c;
        }

        __host__ __device__ const Cylinder_Y &getCylinder_Y() const
        {
            assert(type == CYLINDER_Y);
            return cylinder_y;
        }

        void setRectangleParabolic(const Rectangle_Parabolic &r)
        {
            assert(type == UNKNOWN_TYPE);
            type = RECTANGLE_PARABOLIC;
            rectangle_parabolic = r;
        }

        __host__ __device__ const Rectangle_Parabolic &getRectangleParabolic() const
        {
            assert(type == RECTANGLE_PARABOLIC);
            return rectangle_parabolic;
        }

        void setTriangle_Flat(const Triangle_Flat &t)
        {
            assert(type == UNKNOWN_TYPE);
            type = TRIANGLE_FLAT;
            triangle_flat = t;
        }

        __host__ __device__ const Triangle_Flat &getTriangle_Flat() const
        {
            assert(type == TRIANGLE_FLAT);
            return triangle_flat;
        }

        void setQuadrilateral_Flat(const Quadrilateral_Flat &q)
        {
            assert(type == UNKNOWN_TYPE);
            type = QUADRILATERAL_FLAT;
            quadrilateral_flat = q;
        }

        __host__ __device__ const Quadrilateral_Flat &getQuadrilateral_Flat() const
        {
            assert(type == QUADRILATERAL_FLAT);
            return quadrilateral_flat;
        }

        void setCircle_Flat(const Circle_Flat &c)
        {
            assert(type == UNKNOWN_TYPE);
            type = CIRCLE_FLAT;
            circle_flat = c;
        }

        __host__ __device__ const Circle_Flat &getCircle_Flat() const
        {
            assert(type == CIRCLE_FLAT);
            return circle_flat;
        }

        void setHexagon_Flat(const Hexagon_Flat &h)
        {
            assert(type == UNKNOWN_TYPE);
            type = HEXAGON_FLAT;
            hexagon_flat = h;
        }

        __host__ __device__ const Hexagon_Flat &getHexagon_Flat() const
        {
            assert(type == HEXAGON_FLAT);
            return hexagon_flat;
        }

        void setAnnulus_Flat(const Annulus_Flat &anf)
        {
            assert(type == UNKNOWN_TYPE);
            type = ANNULUS_FLAT;
            annulus_flat = anf;
        }

        __host__ __device__ const Annulus_Flat &getAnnulus_Flat() const
        {
            assert(type == ANNULUS_FLAT);
            return annulus_flat;
        }

        void setCircle_Parabolic(const Circle_Parabolic &circp)
        {
            assert(type == UNKNOWN_TYPE);
            type = CIRCLE_PARABOLIC;
            circle_parabolic = circp;
        }

        __host__ __device__ const Circle_Parabolic &getCircle_Parabolic() const
        {
            assert(type == CIRCLE_PARABOLIC);
            return circle_parabolic;
        }

        void setHexagon_Parabolic(const Hexagon_Parabolic &hexp)
        {
            assert(type == UNKNOWN_TYPE);
            type = HEXAGON_PARABOLIC;
            hexagon_parabolic = hexp;
        }

        __host__ __device__ const Hexagon_Parabolic &getHexagon_Parabolic() const
        {
            assert(type == HEXAGON_PARABOLIC);
            return hexagon_parabolic;
        }

        void setTriangle_Parabolic(const Triangle_Parabolic &tp)
        {
            assert(type == UNKNOWN_TYPE);
            type = TRIANGLE_PARABOLIC;
            triangle_parabolic = tp;
        }

        __host__ __device__ const Triangle_Parabolic &getTriangle_Parabolic() const
        {
            assert(type == TRIANGLE_PARABOLIC);
            return triangle_parabolic;
        }

        void setAnnulus_Parabolic(const Annulus_Parabolic &ap)
        {
            assert(type == UNKNOWN_TYPE);
            type = ANNULUS_PARABOLIC;
            annulus_parabolic = ap;
        }

        __host__ __device__ const Annulus_Parabolic &getAnnulus_Parabolic() const
        {
            assert(type == ANNULUS_PARABOLIC);
            return annulus_parabolic;
        }

        void setQuadrilateral_Parabolic(const Quadrilateral_Parabolic &qp)
        {
            assert(type == UNKNOWN_TYPE);
            type = QUADRILATERAL_PARABOLIC;
            quadrilateral_parabolic = qp;
        }

        __host__ __device__ const Quadrilateral_Parabolic &getQuadrilateral_Parabolic() const
        {
            assert(type == QUADRILATERAL_PARABOLIC);
            return quadrilateral_parabolic;
        }

        void setRectangle_Spherical(const Rectangle_Spherical &rs)
        {
            assert(type == UNKNOWN_TYPE);
            type = RECTANGLE_SPHERICAL;
            rectangle_spherical = rs;
        }

        __host__ __device__ const Rectangle_Spherical &getRectangle_Spherical() const
        {
            assert(type == RECTANGLE_SPHERICAL);
            return rectangle_spherical;
        }

        void setCircle_Spherical(const Circle_Spherical &cs)
        {
            assert(type == UNKNOWN_TYPE);
            type = CIRCLE_SPHERICAL;
            circle_spherical = cs;
        }

        __host__ __device__ const Circle_Spherical &getCircle_Spherical() const
        {
            assert(type == CIRCLE_SPHERICAL);
            return circle_spherical;
        }

        void setHexagon_Spherical(const Hexagon_Spherical &hs)
        {
            assert(type == UNKNOWN_TYPE);
            type = HEXAGON_SPHERICAL;
            hexagon_spherical = hs;
        }

        __host__ __device__ const Hexagon_Spherical &getHexagon_Spherical() const
        {
            assert(type == HEXAGON_SPHERICAL);
            return hexagon_spherical;
        }

        void setAnnulus_Spherical(const Annulus_Spherical &as)
        {
            assert(type == UNKNOWN_TYPE);
            type = ANNULUS_SPHERICAL;
            annulus_spherical = as;
        }

        __host__ __device__ const Annulus_Spherical &getAnnulus_Spherical() const
        {
            assert(type == ANNULUS_SPHERICAL);
            return annulus_spherical;
        }

        void setTriangle_Spherical(const Triangle_Spherical &ts)
        {
            assert(type == UNKNOWN_TYPE);
            type = TRIANGLE_SPHERICAL;
            triangle_spherical = ts;
        }

        __host__ __device__ const Triangle_Spherical &getTriangle_Spherical() const
        {
            assert(type == TRIANGLE_SPHERICAL);
            return triangle_spherical;
        }

        void setQuadrilateral_Spherical(const Quadrilateral_Spherical &qs)
        {
            assert(type == UNKNOWN_TYPE);
            type = QUADRILATERAL_SPHERICAL;
            quadrilateral_spherical = qs;
        }

        __host__ __device__ const Quadrilateral_Spherical &getQuadrilateral_Spherical() const
        {
            assert(type == QUADRILATERAL_SPHERICAL);
            return quadrilateral_spherical;
        }

        Type type = UNKNOWN_TYPE;

        int32_t id = OptixCSP::kElementIdUnassigned;

    private:
        union
        {
            Parallelogram parallelogram;
            Cylinder_Y cylinder_y;
            Rectangle_Parabolic rectangle_parabolic;
            Rectangle_Flat rectangle_flat;
            Triangle_Flat triangle_flat;
            Quadrilateral_Flat quadrilateral_flat;
            Circle_Flat circle_flat;
            Hexagon_Flat hexagon_flat;
            Annulus_Flat annulus_flat;
            Circle_Parabolic circle_parabolic;
            Hexagon_Parabolic hexagon_parabolic;
            Triangle_Parabolic triangle_parabolic;
            Annulus_Parabolic annulus_parabolic;
            Quadrilateral_Parabolic quadrilateral_parabolic;
            Rectangle_Spherical rectangle_spherical;
            Circle_Spherical circle_spherical;
            Hexagon_Spherical hexagon_spherical;
            Annulus_Spherical annulus_spherical;
            Triangle_Spherical triangle_spherical;
            Quadrilateral_Spherical quadrilateral_spherical;
        };
    };
}
