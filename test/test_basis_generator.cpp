//
// Created by Jagan on 06/08/2023.
//
#include "BasisGenerator.hpp"

#include <gtest/gtest.h>

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
  State validState{{{2, 3, 5}, {1, 2, 4}, {3, 2, 1}}, {2, 3}};
  EXPECT_TRUE(basis_generator.check_state(validState));

  // Test case 2: Invalid state due to mismatched qudit length
  State invalidState1{{{1, 2, 3}, {4, 5}, {1, 2, 3}}, {2, 3}};
  EXPECT_FALSE(basis_generator.check_state(invalidState1));

  // Test case 3: Invalid state due to mismatched number of qudits and roots
  State invalidState2{{{1, 2, 3}, {4, 5, 6}}, {2, 3}};
  EXPECT_FALSE(basis_generator.check_state(invalidState2));

  // Test case 4: Invalid state due to failed check_rule
  State invalidState3{{{1, 0, 0}, {4, 5, 6}, {7, 8, 9}}, {0, 5}};
  EXPECT_FALSE(basis_generator.check_state(invalidState3));
}

// Test cases for genState function
TEST(GenStateTest, TestCases) {
  BasisGenerator basis_generator = BasisGenerator();
  // Test case 1: nb_qudits = 3, qudit_len = 2, comb = [1, 2, 3, 4, 5, 6, 7, 8,
  // 9]
  std::vector<int> comb1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  State expectedState1{{{1, 2}, {3, 4}, {5, 6}}, {7, 8, 9}};
  State result1 = basis_generator.gen_state(comb1, 3, 2);
  EXPECT_EQ(result1.qudits, expectedState1.qudits);
  EXPECT_EQ(result1.roots, expectedState1.roots);

  // Test case 2: nb_qudits = 2, qudit_len = 3, comb = [1, 2, 3, 4, 5, 6, 7]
  std::vector<int> comb2 = {1, 2, 3, 4, 5, 6, 7};
  State expectedState2{{{1, 2, 3}, {4, 5, 6}}, {7}};
  State result2 = basis_generator.gen_state(comb2, 2, 3);
  EXPECT_EQ(result2.qudits, expectedState2.qudits);
  EXPECT_EQ(result2.roots, expectedState2.roots);
}

// Test cases for generateBasis function
TEST(GenerateBasisTest, TestCases) {
  BasisGenerator basis_generator = BasisGenerator();
  // Test case 1: nb_qudits = 3, nb_anyons_per_qudit = 2
  int nb_qudits1 = 1;
  int nb_anyons_per_qudit1 = 2;
  Basis expectedBasis1 = {{{{0}}, {}}, {{{1}}, {}}};
  Basis result1 =
      basis_generator.generate_basis(nb_qudits1, nb_anyons_per_qudit1);
  EXPECT_EQ(result1, expectedBasis1);
}
