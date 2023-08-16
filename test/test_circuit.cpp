//
// Created by Jagan on 05/08/2023.
//
#include "Circuit.hpp"

#include <gtest/gtest.h>
#include <iostream>

class TestCircuit : public ::testing::Test {
protected:
  Circuit circuit;
public:
  void SetUp() override { circuit = Circuit(2, 3); }
};

// Tests factorial of 0.
TEST_F(TestCircuit, CircuitTestInit) {
  // H_2
  Sequence hadSequence2 = {{4, 2}, {5, 2},  {4, -2}, {5, -2}, {4, 2},
                             {5, 4}, {4, -2}, {5, 2},  {4, 2},  {5, -2},
                             {4, 2}, {5, -2}, {4, 4}};

  // CNOT_2->1
  Sequence cnotSequence = {
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, -1}, {3, -1}, {3, -1}, {4, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, -1}, {2, -1}, {2, -1}, {1, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, -1}, {2, -1}, {2, -1}, {1, -1}, {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {1, -1}, {2, -1}, {2, -1}, {1, -1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {2, 1},  {3, 1},
      {1, 1},  {2, 1},  {2, 1},  {1, 1},  {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, 1},  {2, 1},  {2, 1},  {1, 1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {1, -1}, {2, -1}, {2, -1}, {1, -1},
      {3, -1}, {2, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {2, 1},  {3, 1},  {3, 1},  {2, 1},  {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, -1}, {3, -1},
      {3, -1}, {2, -1}, {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1},
      {3, -1}, {4, -1}, {2, -1}, {3, -1}, {3, -1}, {2, -1}, {4, 1},  {3, 1},
      {3, 1},  {4, 1},  {4, 1},  {3, 1},  {3, 1},  {4, 1},  {2, 1},  {3, 1},
      {3, 1},  {2, 1},  {4, -1}, {3, -1}, {3, -1}, {4, -1}, {4, -1}, {3, -1}};

  Eigen::VectorXcd initSequence(13);
  initSequence << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  circuit.initialize(initSequence);
  circuit.braidSequence(hadSequence2);
  circuit.braidSequence(cnotSequence);
  circuit.measure();
  const auto res = circuit.run(1000000);

  for (const auto& [result, count] : res.counts_dict) {
    std::cout << result << ": " << count << "\n";
  }
}
