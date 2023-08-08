//
// Created by Jagan on 05/08/2023.
//
#include "gtest/gtest.h"
#include "circuit.cpp"


class TestCircuit : public ::testing::Test {
protected:
    Circuit *circuit;

    virtual void SetUp() {
        circuit = new Circuit(2, 3);
    }

    virtual void TearDown() {
        delete circuit;
    }
};

// Tests factorial of 0.
TEST(CircuitTest, CircuitTestInit) {
    EXPECT_EQ(1, 1);
}