#include <iostream>
#include "circuit.h"
#include "utils.h"
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
    std::vector<std::tuple<int, int>> braids_history;
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

    void get_sigmas() {
        this->sigmas.reserve(nb_anyons - 1); // Reserve space for sigmas

        for (int index = 1; index < nb_anyons; ++index) {
            Eigen::MatrixXcd sigma = this->operator_generator->generate_braiding_operator(
                    index, nb_qudits, nb_anyons_per_qudit
            );
            this->sigmas.push_back(sigma);
        }
    }

    /**
 * Initializes the circuit in the state input_state.
 *
 * @param input_state A normalized quantum state with the same dimensions as the
 *                    fusion space. For example, for a 1-qudit circuit with 3 anyons,
 *                    input_state must be a 3 dimensional vector with norm 1.
 *
 * @throws std::exception Will be thrown if initialization is attempted after performing
 *                       braiding operations.
 * @throws std::invalid_argument Will be thrown if the input state has the wrong dimension or is not normalized.
 *
 * @return A reference to the same circuit.
 */
    CircuitInterface &initialize(const Eigen::VectorXcd &input_state) override {
        if (nb_braids > 0) {
            throw InitializationException("Initialization should happen before any braiding operation is performed!");
        }

        if (input_state.size() != dim) {
            throw std::invalid_argument("The state has wrong dimension. Should be " + std::to_string(dim));
        }

        double norm = input_state.dot(input_state.conjugate()).real();
        if (Eigen::NumTraits<double>::epsilon() <= std::abs(norm - 1)) {
            throw std::invalid_argument("The input state is not normalized correctly!");
        }

        initial_state = input_state;

        return *this;
    }

    CircuitInterface &braid(int n, int m) override {
        if (measured) {
            throw std::runtime_error("System already measured! Cannot perform further braiding!");
        }

        if (m < 1 || n < 1) {
            throw std::invalid_argument("n, m must be higher than 0!");
        }

        if (m > nb_anyons || n > nb_anyons) {
            throw std::invalid_argument(
                    "The system has only " + std::to_string(nb_anyons) + " anyons! n, m are erroneous!");
        }

        if (std::abs(n - m) != 1) {
            throw std::runtime_error("You can only braid adjacent anyons!");
        }

        if (n < m) {
            unitary = sigmas[n - 1] * unitary;
        } else {
            unitary = sigmas[m - 1].adjoint() * unitary;
        }

        this->braids_history.emplace_back(n, m);

        nb_braids++;

        // Assuming drawer.braid function exists
        // drawer.braid(m, n);

        return *this;

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
