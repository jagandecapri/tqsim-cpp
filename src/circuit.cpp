#include <iostream>
#include <map>
#include <random>
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

        this->generate_basis();
        dim = basis.size();
        initial_state = Eigen::VectorXcd::Zero(dim);
        initial_state(0) = 1.0;

        this->get_sigmas();
//        for (int i = 0; i < sigmas.size(); ++i) {
//            std::cout << "Sigmas: " << i << " " << sigmas[i] << std::endl;
//        }
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

    std::vector<std::tuple<int, int>> get_braids_history() override {
        return braids_history;
    }

    vector<Eigen::MatrixXcd> get_braiding_operators() override {
        return sigmas;
    }

    Eigen::MatrixXcd get_unitary() override {
        return unitary;
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
        //std::cout << "Input state: " << input_state << std::endl;
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

//        std::cout << "Unitary before modificaition: " << unitary << std::endl;

        if (n < m) {
            unitary = sigmas[n - 1] * unitary;
        } else {
            unitary = sigmas[m - 1].adjoint() * unitary;
        }

        this->braids_history.emplace_back(n, m);

        nb_braids++;

//        std::cout << "n = " << n << " m = " << m << std::endl;
//        std::cout << "nb_braids: " << nb_braids << std::endl;
//        if (n < m) {
//            std::cout << "n<m sigmas: " << sigmas[n - 1] << std::endl;
//        } else {
//            std::cout << "n >= m sigmas: " << sigmas[m - 1].adjoint() << std::endl;
//        }
//        std::cout << "unitary: " << unitary << std::endl;

        // Assuming drawer.braid function exists
        // drawer.braid(m, n);

        return *this;

    }

/**
     * Takes a sequence of [sigma operator, power], and applies the successive operators to the 'power'.
     * The first operator in the sequence is the first to be applied.
     *
     * @param braid A sequence of pairs of integers representing a braiding operator
     *              and an exponent.
     *
     * @throws std::invalid_argument Is thrown if one of the operators' indices is not an integer
     *                               greater or equal to 1, or is an incorrect index.
     *
     * @return A reference to the same circuit.
     */
    Circuit &braid_sequence(Sequence &braid) override {
        for (const auto &step: braid) {
            if (step.size() != 2 ||
                !std::all_of(step.begin(), step.end(), [](int val) { return std::is_integral<int>::value; })) {
                throw std::invalid_argument("Indices and powers must be integers!");
            }

            int ind = step[0];
            int power = step[1];

            if (ind < 1 || ind >= nb_anyons) {
                throw std::invalid_argument("Invalid braiding operator index!");
            }

            // Computing m and n
            int m = 0, n = 0;
            if (power > 0) {
                n = ind;
                m = ind + 1;
            } else if (power < 0) {
                m = ind;
                n = ind + 1;
            } else {  // if power=0, do nothing (identity)
                continue;
            }

            for (int _ = 0; _ < std::abs(power); ++_) {
                this->braid(n, m); // Assuming the braid function exists
            }
        }
//        std::cout << "Final unitary: " << unitary << std::endl;
        return *this;
    }

    Circuit &measure() override {
        if (this->measured) {
            throw std::runtime_error("Cannot carry the measurements twice!");
        }
        this->measured = true;
        // Assuming drawer.measure function exists
        // drawer.measure();
        return *this;
    }

    // Constructor and other methods

    /**
     * Computes and returns the current state vector of the circuit.
     *
     * @return Eigen::VectorXcd The state vector of the circuit.
     */
    Eigen::VectorXcd statevector() override {
//        std::cout << "Statevector computation";
//        std::cout << "Unitary " << unitary << std::endl;
//        std::cout << "Initial state: " << initial_state << std::endl;
//        std::cout << "Final state: " << unitary * initial_state << std::endl;
        return unitary * initial_state;
    }

    /**
* Simulates the quantum circuit for a specified number of shots and returns the measurement results.
*
* @param shots Number of times the circuit is simulated.
* @return std::map<int, int> Contains the number of measurements for each measured state.
* @throws std::runtime_error Is thrown if the circuit is run without a measurement.
*/
    std::map<int, int> run(int shots) override {
        if (!measured) {
            throw std::runtime_error("The system was not measured!");
        }

        Eigen::VectorXcd statevector = this->statevector();
//        std::cout << "Statevector: " << statevector << std::endl;
        Eigen::VectorXd probs = (statevector.cwiseProduct(statevector.conjugate())).real();
//        std::cout << "Probs: " << probs << std::endl;
        std::discrete_distribution<int> distribution(probs.data(), probs.data() + probs.size());

        // Create a random number generator
        std::random_device rd;
        std::mt19937 generator(rd());

        std::vector<int> memory(shots);
        for (int i = 0; i < shots; ++i) {
            memory[i] = distribution(generator);
        }

        std::map<int, int> counts_dict;
        for (int value: memory) {
            counts_dict[value]++;
        }

        return counts_dict;
    }
};

// int main(int argc, char const *argv[])
// {
//     /* code */
//     Circuit circuit = Circuit(1);
//     circuit.print();
//     return 0;
// }
