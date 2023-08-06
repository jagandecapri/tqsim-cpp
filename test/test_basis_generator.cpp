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

// Test cases for checkState function
TEST(CheckStateTest, TestCases) {
    BasisGenerator basis_generator = BasisGenerator();
    // Test case 1: Valid state
    Basis validState{
            {{2, 3, 5}, {1, 2, 4}, {3, 2, 1}},
            {2,         3}
    };
    EXPECT_TRUE(basis_generator.check_state(validState));

    // Test case 2: Invalid state due to mismatched qudit length
    Basis invalidState1{
            {{1, 2, 3}, {4, 5}, {1, 2, 3}},
            {2,         3}
    };
    EXPECT_FALSE(basis_generator.check_state(invalidState1));

    // Test case 3: Invalid state due to mismatched number of qudits and roots
    Basis invalidState2{
            {{1, 2, 3}, {4, 5, 6}},
            {2,         3}
    };
    EXPECT_FALSE(basis_generator.check_state(invalidState2));

    // Test case 4: Invalid state due to failed check_rule
    Basis invalidState3{
            {{1, 0, 0}, {4, 5, 6}, {7, 8, 9}},
            {0,         5}
    };
    EXPECT_FALSE(basis_generator.check_state(invalidState3));
}