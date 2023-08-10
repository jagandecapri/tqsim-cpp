//
// Created by Jagan on 06/08/2023.
//
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "container.h"
#include "operator_generator.cpp"
#include "basis_generator.cpp"


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

// Test case for the sigma function
TEST(SigmaFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1
    int index_case1 = 1;
    std::vector<int> state_f_case1 = {1, 0};
    std::vector<int> state_i_case1 = {1, 0};

    std::complex<double> result_case1 = operator_generator.sigma(index_case1, state_f_case1, state_i_case1);

    // You need to define the expected result for the test case based on your implementation.
    // For example:
    std::complex<double> expected_case1(-0.30901699437494734, 0.9510565162951536);
    EXPECT_EQ(result_case1, expected_case1);

    // Test case 2: Add more test cases as needed.
}

// Test case for the L function
TEST(LFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1
    int k_case1 = 0;
    int h_case1 = 1;
    int i_case1 = 0;
    int icase1 = 0;
    std::vector<int> jj_case1 = {1, 0};
    std::vector<int> jjcase1 = {1, 0};

    std::complex<double> result_case1 = operator_generator.L(k_case1, h_case1, i_case1, icase1, jj_case1, jj_case1);

    // You need to define the expected result for the test case based on your implementation.
    // For example:
    std::complex<double> expected_case1(-0.5, 0.3632712640026805);
    EXPECT_EQ(result_case1, expected_case1);

    // Test case 2: Add more test cases as needed.
}

// Test case for the S function
TEST(SFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1
    int jm_case1 = 0;
    int jmo_case1 = 1;
    int jmoo_case1 = 0;
    int jmocase1 = 0;
    int h_case1 = 0;
    int i_case1 = 1;
    int icase1 = 0;
    std::vector<int> jj_case1 = {0, 1};
    std::vector<int> jjcase1 = {1, 0};

    std::complex<double> result_case1 = operator_generator.S(jm_case1, jmo_case1, jmoo_case1, jmocase1, h_case1,
                                                             icase1, i_case1, jjcase1, jj_case1);

    // You need to define the expected result for the test case based on your implementation.
    // For example:
    std::complex<double> expected_case1(0, 0);
    EXPECT_EQ(result_case1, expected_case1);

    // Test case 2: Add more test cases as needed.
}

TEST(GenSigmaFunctionTest, TestCases) {
    OperatorGenerator operator_generator = OperatorGenerator();
    // Test case 1
    int index_case1 = 3;

    State state_i_case_1 = {
            .qudits =  {{0, 1},
                        {0, 1}},
            .roots =  {0}
    };

    State state_f_case_1 = {
            .qudits =  {{1, 0},
                        {1, 0}},
            .roots =  {0}
    };

    std::complex<double> result_case1 = operator_generator.gen_sigma(index_case1, state_i_case_1, state_f_case_1);
    std::complex<double> expected_case1(-0, 0);
    EXPECT_EQ(result_case1, expected_case1);
}

TEST(GenerateBraidingOperator, TestCases) {
    BasisGeneratorInterface *basis_generator = new BasisGenerator();
    OperatorGenerator operator_generator = OperatorGenerator(basis_generator);

    int index_case1 = 1;
    int nb_qudits_case1 = 1;
    int nb_anyons_per_qudit_case1 = 2;

    Eigen::MatrixXcd expected_case1(2, 2);
    expected_case1 << std::complex<double>(-0.8090169943749473, -0.5877852522924732), std::complex<double>(-0, 0),
            std::complex<double>(0, 0), std::complex<double>(-0.30901699437494734, 0.9510565162951536);

    Eigen::MatrixXcd result_case1 = operator_generator.generate_braiding_operator(index_case1, nb_qudits_case1,
                                                                                  nb_anyons_per_qudit_case1);
    EXPECT_EQ(result_case1, expected_case1);
}

TEST(GenerateBraidingOperator1, TestCases) {
    BasisGeneratorInterface *basis_generator = new BasisGenerator();
    OperatorGenerator operator_generator = OperatorGenerator(basis_generator);

    int index_case1 = 3;
    int nb_qudits_case1 = 2;
    int nb_anyons_per_qudit_case1 = 3;

//    Eigen::MatrixXcd expected_case1(2, 2);
//    expected_case1 << std::complex<double>(-0.8090169943749473, -0.5877852522924732), std::complex<double>(-0, 0),
//            std::complex<double>(0, 0), std::complex<double>(-0.30901699437494734, 0.9510565162951536);

    Eigen::MatrixXcd result_case1 = operator_generator.generate_braiding_operator(index_case1, nb_qudits_case1,
                                                                                  nb_anyons_per_qudit_case1);
    std::cout << result_case1 << std::endl;

//    EXPECT_EQ(result_case1, expected_case1);
}