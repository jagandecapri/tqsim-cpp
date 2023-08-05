#include <iostream>
#include "circuit.h"

using namespace std;

Circuit::Circuit(int nb_qudits, int nb_anyons_per_qudit)
{
    this->nb_qudits = nb_qudits;
    this->nb_anyons_per_qudit = nb_anyons_per_qudit;
}

Circuit::~Circuit()
{
}

void Circuit::print(){
    cout << "nb_qudits = " << nb_qudits <<  " nb_anyons_per_qudit = " << nb_anyons_per_qudit << endl;
};

// int main(int argc, char const *argv[])
// {
//     /* code */
//     Circuit circuit = Circuit(1);
//     circuit.print();
//     return 0;
// }
