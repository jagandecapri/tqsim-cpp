#include <iostream>
#include "circuit.h"
#include "operator_generator.cpp"
#include "basis_generator.cpp"

using namespace std;

class Circuit : public CircuitInterface {
private:
    /* data */
    OperatorGeneratorInterface *operator_generator;
    BasisGeneratorInterface *basis_generator;
    int nb_qudits;
    int nb_anyons_per_qudit;
    int nb_anyons;
    int nb_braids;
    bool measured;
    Basis basis;
    size_t dim;
    Eigen::VectorXcd initial_state;
    Sigma sigmas;
    Eigen::MatrixXcd unitary;
public:
    Circuit(int nb_qudits, int nb_anyons_per_qudit) :
            nb_qudits(nb_qudits), nb_anyons_per_qudit(nb_anyons_per_qudit),
            nb_anyons(nb_qudits * nb_anyons_per_qudit), nb_braids(0),
            measured(false) {
        basis_generator = new BasisGenerator();
        operator_generator = new OperatorGenerator(basis_generator);
        basis = generate_basis();
        dim = basis.size();
        initial_state = Eigen::VectorXcd::Zero(dim);
        initial_state(0) = 1.0;
        sigmas = get_sigmas();
        unitary = Eigen::MatrixXcd::Identity(dim, dim);
    }

    int get_nb_qudits() override {
        return nb_qudits;
    }

    int get_nb_anyons_per_qudit() override {
        return nb_anyons_per_qudit;
    }

    int get_dim() override {
        return dim;
    }

    Basis get_basis() override {
        return basis;
    }

    Sigma get_braiding_operators() override {
        return sigmas;
    }

    Basis generate_basis() {
        basis = basis_generator->generate_basis(nb_qudits, nb_anyons_per_qudit);
        return basis;
    }

    Sigma get_sigmas() {

    }

    void initialize() override {
    }

    void braid(int i, int j) override {

    }

    void braid_sequence(Sequence &sequence) override {

    }

    void measure() override {

    }

    std::string history(std::string) override {

    }

    Eigen::MatrixXd statevector() override {

    }

    void run(int) override {

    }

    void print() {
        cout << "nb_qudits = " << nb_qudits << " nb_anyons_per_qudit = " << nb_anyons_per_qudit << endl;
    }
};

// int main(int argc, char const *argv[])
// {
//     /* code */
//     Circuit circuit = Circuit(1);
//     circuit.print();
//     return 0;
// }
