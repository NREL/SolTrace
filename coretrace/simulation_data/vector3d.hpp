/**
 * @file vector3d.hpp
 * @brief 3D vector and matrix classes
 *
 * Provides Vector3d and Matrix3d classes for 3D geometric operations,
 * coordinate transformations, and linear algebra computations.
 * These classes form the mathematical foundation for ray tracing
 * calculations and geometric transformations in SolTrace.
 *
 * @defgroup math Mathematical Operations
 * @{
 */

#ifndef SOLTRACE_VECTOR3D_H
#define SOLTRACE_VECTOR3D_H

#include <cassert>
#include <vector>

#include "matvec.hpp"

namespace SolTrace::Data {

class Vector3d
{
public:
    /**
     * @brief Default constructor - initializes to zero vector
     */
    Vector3d();

    /**
     * @brief Constructor from array
     * @param data Array of 3 doubles [x, y, z]
     */
    Vector3d(const double data[3]);

    /**
     * @brief Constructor from components
     * @param x X component
     * @param y Y component
     * @param z Z component
     */
    Vector3d(double x, double y, double z);
    ~Vector3d();

    /**
     * @brief Set all components to zero
     */
    void zero();

    /**
     * @brief Set vector components
     * @param x X component
     * @param y Y component
     * @param z Z component
     */
    void set_values(double x, double y, double z);

    /**
     * @brief Multiply vector by scalar
     * @param alpha Scalar multiplier
     */
    void scalar_mult(double alpha);

    /**
     * @brief Make the vector a unit vector
     */
    void make_unit();

    /**
     * @brief Compute Euclidean vector norm
     * @return Euclidean norm
     */
    double norm() const;

    /**
     * @brief Const access operator
     * @param idx Component index (0=x, 1=y, 2=z)
     * @return Const reference to component
     */
    const double &operator[](int idx) const;

    /**
     * @brief Access operator
     * @param idx Component index (0=x, 1=y, 2=z)
     * @return Reference to component
     */
    double &operator[](int idx);

    /**
     * @brief Stream output operator
     * @param os Output stream
     * @param x Vector to output
     * @return Reference to output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const Vector3d &x);

    double data[3];

private:
};

class Matrix3d
{
public:
    /**
     * @brief Default constructor - initializes to zero matrix
     */
    Matrix3d();
    ~Matrix3d();

    /**
     * @brief Set matrix element value
     * @param i Row index (0-2)
     * @param j Column index (0-2)
     * @param val Value to set
     */
    void set_value(int i, int j, double val);

    /**
     * @brief Get matrix element value
     * @param i Row index (0-2)
     * @param j Column index (0-2)
     * @return Matrix element value
     */
    double get_value(int i, int j) const;

    /**
     * @brief Set all matrix elements to zero
     */
    void zero();

    /**
     * @brief Set matrix to identity (diagonal = 1, off-diagonal = 0)
     */
    void identity();

    // double data[9];
    double data[3][3];

    friend std::ostream &operator<<(std::ostream &os, const Matrix3d &A);

private:
};

inline void vector_copy(double data[3], const Vector3d &x)
{
    CopyVec3(data, x.data);
    return;
}
inline void vector_copy(std::vector<double> &dest, const Vector3d &x)
{
    CopyVec3(dest, x.data);
    return;
}

inline void matrix_copy(double data[3][3], const Matrix3d &A)
{
    CopyMat3(data, A.data);
    return;
}

// Compute y = A*x placing the result in y
void matrix_vector_product(const Matrix3d &A, const Vector3d &x, Vector3d &y);
// Compute C = A * B placing result in C
void matrix_matrix_product(const Matrix3d &A, const Matrix3d &B, Matrix3d &C);

// Compute z = a*x + b*y (result stored in z)
void vector_add(double a, const Vector3d &x,
                double b, const Vector3d &y,
                Vector3d &z);
/**
 * @brief Compute linear combination y = a*x + b*y (result stored in y)
 * @param a Scalar multiplier for vector x
 * @param x First vector
 * @param b Scalar multiplier for vector y
 * @param y Second vector (modified in-place with result)
 */
void vector_add(double a, const Vector3d &x,
                double b, Vector3d &y);

/**
 * @brief Compute element-wise maximum of two vectors
 * @param x First vector
 * @param y Second vector
 * @param max Output vector containing element-wise maximum
 */
void vector_max(const Vector3d &x, const Vector3d &y, Vector3d &max);

/**
 * @brief Compute element-wise minimum of two vectors
 * @param x First vector
 * @param y Second vector
 * @param min Output vector containing element-wise minimum
 */
void vector_min(const Vector3d &x, const Vector3d &y, Vector3d &min);

/**
 * @brief Compute standard Euclidean dot product
 * @param x First vector
 * @param y Second vector
 * @return Dot product x·y
 */
double dot_product(const Vector3d &x, const Vector3d &y);

/**
 * @brief Compute standard Euclidean cross product
 * @param u First vector
 * @param v Second vector
 * @param w Result vector
 */
void cross_product(const Vector3d &u, const Vector3d &v, Vector3d &w);

double error(const Vector3d &u, const Vector3d &v);
double error_inf(const Vector3d &u, const Vector3d &v);

/**
 * @brief Compute Euclidean norm (length) of vector
 * @param x Input vector
 * @return ||x||₂ (L2 norm)
 */
double vector_norm(const Vector3d &x);

/**
 * @brief Normalize vector to unit length (modifies vector in-place)
 * @param x Vector to normalize
 */
void make_unit_vector(Vector3d &x);
// void transform_to_local(Vector3d &pos_ref,
//                         Vector3d &cos_ref,
//                         Vector3d &origin,
//                         Matrix3d &ref_to_local,
//                         Vector3d &pos_local,
//                         Vector3d &cos_local);
// void transform_to_reference(Vector3d &pos_local,
//                             Vector3d &cos_local,
//                             Vector3d &origin,
//                             Matrix3d &local_to_ref,
//                             Vector3d &pos_ref,
//                             Vector3d &cos_ref);

/**
 * @brief Compute transformation matrices from Euler angles
 * @param euler Euler angles vector [rotation_x, rotation_y, rotation_z] in radians
 * @param ref_to_local Output matrix for reference-to-local transformation
 * @param local_to_ref Output matrix for local-to-reference transformation
 */
void compute_transform_matrices(Vector3d &euler,
                                Matrix3d &ref_to_local,
                                Matrix3d &local_to_ref);

} // namespace SolTrace::Data

/**
 * @}
 */

#endif
