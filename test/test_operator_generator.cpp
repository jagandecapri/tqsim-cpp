//
// Created by Jagan on 06/08/2023.
//
#include <iostream>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <operator_generator.h>

// Test case for the F function
TEST(FunctionTest, FTest) {
    // Test case 1: a1 + a2 + a3 + outcome = 4
    int a1_case1 = 1, a2_case1 = 1, a3_case1 = 1, outcome_case1 = 1;
    Eigen::Matrix<double, 2, 2> result_case1 = F(a1_case1, a2_case1, a3_case1, outcome_case1);
    Eigen::Matrix<double, 2, 2> expected_case1;
    const double inv_phi = (std::sqrt(5) - 1) / 2;
    expected_case1 << inv_phi, std::sqrt(inv_phi),
            std::sqrt(inv_phi), -inv_phi;
    ASSERT_TRUE(result_case1.isApprox(expected_case1, 1e-6));

    // Test case 2: a1 + a2 + a3 + outcome = 3
    int a1_case2 = 0, a2_case2 = 1, a3_case2 = 2, outcome_case2 = 0;
    Eigen::Matrix<double, 2, 2> result_case2 = F(a1_case2, a2_case2, a3_case2, outcome_case2);
    Eigen::Matrix<double, 2, 2> expected_case2;
    expected_case2 << 0, 0,
            0, 1;
    ASSERT_TRUE(result_case2.isApprox(expected_case2, 1e-6));

    // Add more test cases as needed...
}
