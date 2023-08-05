//
// Created by Jagan on 05/08/2023.
//
#include "gtest/gtest.h"
#include "circuit.h"


class TestCircuit : public ::testing::Test {
protected:
    virtual void SetUp()
    {
        circuit = new Circuit(1);
    }

    virtual void TearDown() {
        delete circuit;
    }

    Circuit * circuit;
};

// Tests factorial of 0.
TEST(CircuitTest, CircuitTestInit) {
    EXPECT_EQ(1, 1);
}