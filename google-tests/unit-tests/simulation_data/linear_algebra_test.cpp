#include <gtest/gtest.h>

#include <iomanip>
#include <cmath>
#include <sstream>

#include <vector3d.hpp>

#include "common.hpp"

TEST(LinearAlgebra, VectorBasics)
{
    const double TOL = 1e-12;

    Vector3d x(1.0, 1.0, 3.0);
    Vector3d y(-1.0, 2.0, -1.0);
    Vector3d z;

    EXPECT_EQ(z[0], 0.0);
    EXPECT_EQ(z[1], 0.0);
    EXPECT_EQ(z[2], 0.0);

    vector_add(1.0, x, 1.0, y, z);
    EXPECT_NEAR(z[0], 0.0, TOL);
    EXPECT_NEAR(z[1], 3.0, TOL);
    EXPECT_NEAR(z[2], 2.0, TOL);

    vector_add(0.5, x, -2.0, y, z);
    EXPECT_NEAR(z[0], 2.5, TOL);
    EXPECT_NEAR(z[1], -3.5, TOL);
    EXPECT_NEAR(z[2], 3.5, TOL);

    vector_max(x, y, z);
    EXPECT_EQ(z[0], 1.0);
    EXPECT_EQ(z[1], 2.0);
    EXPECT_EQ(z[2], 3.0);

    vector_min(x, y, z);
    EXPECT_EQ(z[0], -1.0);
    EXPECT_EQ(z[1], 1.0);
    EXPECT_EQ(z[2], -1.0);

    z.set_values(1.0, 2.0, 3.0);
    z.scalar_mult(-0.5);
    EXPECT_NEAR(z[0], -0.5, TOL);
    EXPECT_NEAR(z[1], -1.0, TOL);
    EXPECT_NEAR(z[2], -1.5, TOL);

    double dp = dot_product(x, y);
    EXPECT_NEAR(dp, -2.0, TOL);

    double mag = vector_norm(x);
    EXPECT_NEAR(mag, sqrt(11.0), TOL);

    make_unit_vector(x);
    EXPECT_NEAR(x[0], 1.0 / mag, TOL);
    EXPECT_NEAR(x[1], 1.0 / mag, TOL);
    EXPECT_NEAR(x[2], 3.0 / mag, TOL);

    mag = vector_norm(y);
    EXPECT_NEAR(mag, sqrt(6.0), TOL);

    make_unit_vector(y);
    EXPECT_NEAR(y[0], -1.0 / mag, TOL);
    EXPECT_NEAR(y[1], 2.0 / mag, TOL);
    EXPECT_NEAR(y[2], -1.0 / mag, TOL);

    z.zero();
    EXPECT_EQ(z[0], 0.0);
    EXPECT_EQ(z[1], 0.0);
    EXPECT_EQ(z[2], 0.0);

    z[0] = 1.0;
    z[1] = 2.0;
    z[2] = 3.0;
    EXPECT_EQ(z[0], 1.0);
    EXPECT_EQ(z[1], 2.0);
    EXPECT_EQ(z[2], 3.0);

    Vector3d u(z.data);
    EXPECT_EQ(z[0], u[0]);
    EXPECT_EQ(z[1], u[1]);
    EXPECT_EQ(z[2], u[2]);

    Vector3d ihat(1.0, 0.0, 0.0);
    Vector3d jhat(0.0, 1.0, 0.0);
    Vector3d khat(0.0, 0.0, 1.0);
    Vector3d result;
    cross_product(ihat, jhat, result);
    EXPECT_NEAR(result[0], 0.0, TOL);
    EXPECT_NEAR(result[1], 0.0, TOL);
    EXPECT_NEAR(result[2], 1.0, TOL);
    result.zero();

    cross_product(ihat, khat, result);
    EXPECT_NEAR(result[0], 0.0, TOL);
    EXPECT_NEAR(result[1], -1.0, TOL);
    EXPECT_NEAR(result[2], 0.0, TOL);
    result.zero();

    cross_product(jhat, khat, result);
    EXPECT_NEAR(result[0], 1.0, TOL);
    EXPECT_NEAR(result[1], 0.0, TOL);
    EXPECT_NEAR(result[2], 0.0, TOL);
    result.zero();

    double l2 = error(ihat, jhat);
    EXPECT_NEAR(l2, sqrt(2.0), TOL);
    double linf = error_inf(ihat, jhat);
    EXPECT_NEAR(linf, 1.0, TOL);
}

TEST(LinearAlgebra, MatrixVectorProduct)
{
    Vector3d x(0.5, 1.0, -2.0);
    Vector3d y;

    Matrix3d A;
    A.set_value(0, 0, 1.0);
    A.set_value(1, 1, 1.0);
    A.set_value(2, 2, 1.0);

    matrix_vector_product(A, x, y);
    EXPECT_TRUE(is_identical(x, y));

    A.set_value(0, 1, 2.0);
    A.set_value(1, 2, 3.0);
    matrix_vector_product(A, x, y);
    EXPECT_NEAR(y[0], 2.5, 1e-12);
    EXPECT_NEAR(y[1], -5.0, 1e-12);
    EXPECT_NEAR(y[2], -2.0, 1e-12);
}

TEST(LinearAlgebra, MatrixMatrixProduct)
{
    Matrix3d A;
    A.set_value(0, 1, 1.0);
    A.set_value(0, 2, 2.0);
    A.set_value(1, 0, -1.0);
    A.set_value(1, 1, 1.0);
    A.set_value(2, 0, 2.0);
    A.set_value(2, 2, 1.0);

    Matrix3d B;
    B.set_value(0, 0, -1.0);
    B.set_value(0, 1, 1.0);
    B.set_value(1, 1, 1.0);
    B.set_value(1, 2, -2.0);
    B.set_value(2, 2, 1.0);

    Matrix3d Ctrue;
    Ctrue.set_value(0, 1, 1.0);
    Ctrue.set_value(1, 0, 1.0);
    Ctrue.set_value(1, 2, -2.0);
    Ctrue.set_value(2, 0, -2.0);
    Ctrue.set_value(2, 1, 2.0);
    Ctrue.set_value(2, 2, 1.0);

    Matrix3d C;
    C.zero();

    matrix_matrix_product(A, B, C);
    EXPECT_TRUE(is_identical(C, Ctrue));
}

TEST(LinearAlgebra, CoordinateTransforms)
{
    // TODO: Implement tests for compute_transform_matrices,
    // transform_to_local, and transform_to_reference functions
    Matrix3d A;
    Matrix3d B;

    EXPECT_TRUE(is_identical(A, B));
}

TEST(LinearAlgebra, Vector3dOutputOperator)
{
    Vector3d v1(1.0, 2.0, 3.0);
    Vector3d v2(-0.5, 0.0, 10.5);
    Vector3d v3(0.0, 0.0, 0.0);

    std::ostringstream oss1;
    oss1 << v1;
    EXPECT_EQ(oss1.str(), "[1, 2, 3]");

    std::ostringstream oss2;
    oss2 << v2;
    EXPECT_EQ(oss2.str(), "[-0.5, 0, 10.5]");

    std::ostringstream oss3;
    oss3 << v3;
    EXPECT_EQ(oss3.str(), "[0, 0, 0]");

    // Test with different precision
    Vector3d v4(1.23456789, -2.34567891, 3.45678912);
    std::ostringstream oss4;
    oss4 << std::fixed << std::setprecision(3) << v4;
    EXPECT_EQ(oss4.str(), "[1.235, -2.346, 3.457]");
}

TEST(LinearAlgebra, Matrix3dOutputOperator)
{
    Matrix3d A;
    A.zero();
    A.set_value(0, 0, 1.0);
    A.set_value(0, 1, 2.0);
    A.set_value(0, 2, 3.0);
    A.set_value(1, 0, 4.0);
    A.set_value(1, 1, 5.0);
    A.set_value(1, 2, 6.0);
    A.set_value(2, 0, 7.0);
    A.set_value(2, 1, 8.0);
    A.set_value(2, 2, 9.0);

    std::ostringstream oss;
    oss << A;
    EXPECT_EQ(oss.str(), "[1, 2, 3; 4, 5, 6; 7, 8, 9]");

    // Test identity matrix
    Matrix3d I;
    I.identity();
    std::ostringstream oss_identity;
    oss_identity << I;
    EXPECT_EQ(oss_identity.str(), "[1, 0, 0; 0, 1, 0; 0, 0, 1]");

    // Test zero matrix
    Matrix3d Z;
    Z.zero();
    std::ostringstream oss_zero;
    oss_zero << Z;
    EXPECT_EQ(oss_zero.str(), "[0, 0, 0; 0, 0, 0; 0, 0, 0]");

    // Test with negative values and different precision
    Matrix3d B;
    B.set_value(0, 0, -1.5);
    B.set_value(0, 1, 2.25);
    B.set_value(0, 2, -3.75);
    B.set_value(1, 0, 0.0);
    B.set_value(1, 1, -0.5);
    B.set_value(1, 2, 1.0);
    B.set_value(2, 0, 10.0);
    B.set_value(2, 1, -20.0);
    B.set_value(2, 2, 30.0);

    std::ostringstream oss_negative;
    oss_negative << std::fixed << std::setprecision(2) << B;
    EXPECT_EQ(oss_negative.str(), "[-1.50, 2.25, -3.75; 0.00, -0.50, 1.00; 10.00, -20.00, 30.00]");
}

TEST(LinearAlgebra, Matrix3dGetValue)
{
    Matrix3d A;
    A.zero();
    
    // Test getting values from zero matrix
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_EQ(A.get_value(i, j), 0.0);
        }
    }

    // Set specific values and test retrieval
    A.set_value(0, 0, 1.5);
    A.set_value(0, 1, -2.5);
    A.set_value(0, 2, 3.7);
    A.set_value(1, 0, -4.2);
    A.set_value(1, 1, 5.8);
    A.set_value(1, 2, -6.1);
    A.set_value(2, 0, 7.9);
    A.set_value(2, 1, -8.4);
    A.set_value(2, 2, 9.6);

    EXPECT_DOUBLE_EQ(A.get_value(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(A.get_value(0, 1), -2.5);
    EXPECT_DOUBLE_EQ(A.get_value(0, 2), 3.7);
    EXPECT_DOUBLE_EQ(A.get_value(1, 0), -4.2);
    EXPECT_DOUBLE_EQ(A.get_value(1, 1), 5.8);
    EXPECT_DOUBLE_EQ(A.get_value(1, 2), -6.1);
    EXPECT_DOUBLE_EQ(A.get_value(2, 0), 7.9);
    EXPECT_DOUBLE_EQ(A.get_value(2, 1), -8.4);
    EXPECT_DOUBLE_EQ(A.get_value(2, 2), 9.6);

    // Test identity matrix
    Matrix3d I;
    I.identity();
    EXPECT_DOUBLE_EQ(I.get_value(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(I.get_value(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(I.get_value(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(I.get_value(2, 2), 1.0);

    // Test very small and very large values
    A.set_value(0, 0, 1e-15);
    A.set_value(1, 1, 1e15);
    A.set_value(2, 2, -1e-10);
    
    EXPECT_DOUBLE_EQ(A.get_value(0, 0), 1e-15);
    EXPECT_DOUBLE_EQ(A.get_value(1, 1), 1e15);
    EXPECT_DOUBLE_EQ(A.get_value(2, 2), -1e-10);
}
