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

    virtual std::vector<std::tuple<int, int>> get_braids_history() = 0;

    virtual std::vector<Eigen::MatrixXcd> get_braiding_operators() = 0;

    virtual Eigen::MatrixXcd get_unitary() = 0;

    virtual CircuitInterface &initialize(const Eigen::VectorXcd &) = 0;

    virtual CircuitInterface &braid(int, int) = 0;

    virtual CircuitInterface &braid_sequence(Sequence &) = 0;

    virtual CircuitInterface &measure() = 0;

    virtual Eigen::VectorXcd statevector() = 0;

    virtual Result run(int) = 0;
};

#endif // CIRCUIT_H