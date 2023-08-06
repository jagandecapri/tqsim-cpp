//
// Created by Jagan on 06/08/2023.
//
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <operator_generator.cpp>


/**
 * @brief Custom comparison function for complex numbers.
 *
 * @param a Complex number a.
 * @param b Complex number b.
 * @param tol Tolerance for comparison.
 * @return true if both real and imaginary parts are within the tolerance, false otherwise.
 */
bool complexIsApprox(const std::complex<double> &a, const std::complex<double> &b, double tol = 1e-6) {
    return abs(a.real() - b.real()) < tol && abs(a.imag() - b.imag()) < tol;
}

// Test case for the F function
TEST(FFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1: a1 + a2 + a3 + outcome = 4
    int a1_case1 = 1, a2_case1 = 1, a3_case1 = 1, outcome_case1 = 1;
    Eigen::Matrix<double, 2, 2> result_case1 = operator_generator.F(a1_case1, a2_case1, a3_case1, outcome_case1);
    Eigen::Matrix<double, 2, 2> expected_case1;
    const double inv_phi = (std::sqrt(5) - 1) / 2;
    expected_case1 << inv_phi, std::sqrt(inv_phi),
            std::sqrt(inv_phi), -inv_phi;
    ASSERT_TRUE(result_case1.isApprox(expected_case1, 1e-6));

    // Test case 2: a1 + a2 + a3 + outcome = 3
    int a1_case2 = 0, a2_case2 = 1, a3_case2 = 2, outcome_case2 = 0;
    Eigen::Matrix<double, 2, 2> result_case2 = operator_generator.F(a1_case2, a2_case2, a3_case2, outcome_case2);
    Eigen::Matrix<double, 2, 2> expected_case2;
    expected_case2 << 0, 0,
            0, 1;
    ASSERT_TRUE(result_case2.isApprox(expected_case2, 1e-6));

    // Add more test cases as needed...
}

// Test case for the R function
TEST(RFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1: a1 + a2 = 2
    int a1_case1 = 1, a2_case1 = 1;
    Eigen::Matrix<std::complex<double>, 2, 2> result_case1 = operator_generator.R(a1_case1, a2_case1);
    Eigen::Matrix<std::complex<double>, 2, 2> expected_case1;
    const double pi = std::acos(-1);
    std::complex<double> exp_val1 = std::exp(std::complex<double>(0.0, -4 * pi / 5.0));
    std::complex<double> exp_val2 = std::exp(std::complex<double>(0.0, 3 * pi / 5.0));
    expected_case1 << exp_val1, std::complex<double>(0, 0),
            std::complex<double>(0, 0), exp_val2;
    ASSERT_TRUE(complexIsApprox(result_case1(0, 0), expected_case1(0, 0)) &&
                complexIsApprox(result_case1(0, 1), expected_case1(0, 1)) &&
                complexIsApprox(result_case1(1, 0), expected_case1(1, 0)) &&
                complexIsApprox(result_case1(1, 1), expected_case1(1, 1)));

    // Test case 2: a1 + a2 != 2
    int a1_case2 = 0, a2_case2 = 1;
    Eigen::Matrix<std::complex<double>, 2, 2> result_case2 = operator_generator.R(a1_case2, a2_case2);
    Eigen::Matrix<std::complex<double>, 2, 2> expected_case2;
    expected_case2 << std::complex<double>(1, 0), std::complex<double>(0, 0),
            std::complex<double>(0, 0), std::complex<double>(1, 0);
    ASSERT_TRUE(complexIsApprox(result_case2(0, 0), expected_case2(0, 0)) &&
                complexIsApprox(result_case2(0, 1), expected_case2(0, 1)) &&
                complexIsApprox(result_case2(1, 0), expected_case2(1, 0)) &&
                complexIsApprox(result_case2(1, 1), expected_case2(1, 1)));

    // Add more test cases as needed...
}

// Test case for the B function
TEST(BFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1: a1 + a2 + a3 + outcome = 4
    int a0 = 1, a1 = 1, a2 = 1, outcome = 1;
    Eigen::Matrix<std::complex<double>, 2, 2> result_case1 = operator_generator.B(a0, a1, a2, outcome);
    Eigen::Matrix<std::complex<double>, 2, 2> expected_case1;
    expected_case1 << std::complex<double>(-0.5, 0.36327126400268051), std::complex<double>(-0.24293413587832285,
                                                                                            -0.74767439061061047),
            std::complex<double>(-0.24293413587832285, -0.74767439061061047), std::complex<double>(-0.61803398874989468,
                                                                                                   5.5511151231257827E-17);
    ASSERT_TRUE(complexIsApprox(result_case1(0, 0), expected_case1(0, 0)) &&
                complexIsApprox(result_case1(0, 1), expected_case1(0, 1)) &&
                complexIsApprox(result_case1(1, 0), expected_case1(1, 0)) &&
                complexIsApprox(result_case1(1, 1), expected_case1(1, 1)));
}