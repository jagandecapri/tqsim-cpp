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
    std::vector<Eigen::MatrixXcd> sigmas;
    Eigen::MatrixXcd unitary;
public:
    Circuit(int nb_qudits, int nb_anyons_per_qudit) :
            nb_qudits(nb_qudits), nb_anyons_per_qudit(nb_anyons_per_qudit),
            nb_anyons(nb_qudits * nb_anyons_per_qudit), nb_braids(0),
            measured(false) {
        basis_generator = new BasisGenerator();
        operator_generator = new OperatorGenerator(basis_generator);
        dim = basis.size();
        initial_state = Eigen::VectorXcd::Zero(dim);
        initial_state(0) = 1.0;
        unitary = Eigen::MatrixXcd::Identity(dim, dim);
        this->generate_basis();
        this->get_sigmas();
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

    vector<Eigen::MatrixXcd> get_braiding_operators() override {
        return sigmas;
    }

    void generate_basis() {
        this->basis = basis_generator->generate_basis(nb_qudits, nb_anyons_per_qudit);
    }

    Sigma get_sigmas() {
        std::vector<Eigen::MatrixXcd> sigmas;
        sigmas.reserve(nb_anyons - 1); // Reserve space for sigmas

        for (int index = 1; index < nb_anyons; ++index) {
            Eigen::MatrixXcd sigma = this->operator_generator->generate_braiding_operator(
                    index, nb_qudits, nb_anyons_per_qudit
            );
            sigmas.push_back(sigma);
        }
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
