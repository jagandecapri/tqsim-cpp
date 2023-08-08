#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <Eigen/Dense>
#include "container.h"
#include "operator_generator.h"
#include "basis_generator.h"
// Path: src/circuit.cpp

class CircuitInterface {
public:

    virtual int get_nb_qudits() = 0;

    virtual int get_nb_anyons_per_qudit() = 0;

    virtual int get_dim() = 0;

    virtual Basis get_basis() = 0;

    virtual std::vector<Eigen::MatrixXcd> get_braiding_operators() = 0;

    virtual void initialize() = 0;

    virtual void braid(int, int) = 0;

    virtual void braid_sequence(Sequence &) = 0;

    virtual void measure() = 0;

    virtual std::string history(std::string) = 0;

    virtual Eigen::MatrixXd statevector() = 0;

    virtual void run(int) = 0;

    void print();
};

#endif // CIRCUIT_H