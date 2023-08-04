#ifndef CIRCUIT_H
#define CIRCUIT_H

// Path: src/circuit.cpp

class Circuit{
    private:
        /* data */
        int x;
    public:
        Circuit(int);
        ~Circuit();
        void print();
};

#endif // CIRCUIT_H