#ifndef CIRCUIT_H
#define CIRCUIT_H

// Path: src/circuit.cpp

class Circuit{
    private:
        /* data */
        int nb_qudits;
        int nb_anyons_per_qudit;
    public:
        Circuit(int nb_qudits = 1, int nb_anyons_per_qudit = 3);
        ~Circuit();
        void print();
};

#endif // CIRCUIT_H