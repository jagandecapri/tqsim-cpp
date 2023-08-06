//
// Created by Jagan on 06/08/2023.
//
#include "gtest/gtest.h"
#include "basis_generator.cpp"

// Test cases for checkFibonacciFusionRules function
TEST(CheckRulesTest, TestCases) {
    BasisGenerator basis_generator = BasisGenerator();
    // Test case 1: anyon1 = 2, anyon2 = 3, outcome = 5
    EXPECT_TRUE(basis_generator.check_rule(1, 1, 0));

    // Test case 2: anyon1 = 0, anyon2 = 1, outcome = 1
    EXPECT_TRUE(basis_generator.check_rule(0, 1, 1));

    // Test case 3: anyon1 = 1, anyon2 = 0, outcome = 1
    EXPECT_TRUE(basis_generator.check_rule(1, 0, 1));

    // Test case 4: anyon1 = 0, anyon2 = 0, outcome = 0
    EXPECT_TRUE(basis_generator.check_rule(0, 0, 0));

    // Test case 4: anyon1 = 1, anyon2 = 0, outcome = 0
    EXPECT_FALSE(basis_generator.check_rule(1, 0, 0));
}

TEST(CheckOutcomesTest, TestCases) {
    BasisGenerator basis_generator = BasisGenerator();
    // Test case 1: outcomes = [1, 1, 0, 0, 0]
    std::vector<int> outcomes = {1, 1, 1, 1, 1};
    EXPECT_TRUE(basis_generator.check_outcomes(outcomes));

    // Test case 2: outcomes = [1, 1, 0, 0, 1]
    outcomes = {1, 1, 0, 0, 1};
    EXPECT_FALSE(basis_generator.check_outcomes(outcomes));

    // Test case 3: outcomes = [1, 1, 0, 0, 0, 0]
    outcomes = {1, 1, 0, 0, 0, 0};
    EXPECT_FALSE(basis_generator.check_outcomes(outcomes));
}