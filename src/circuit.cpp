#include <iostream>
#include "circuit.h"

using namespace std;

// class Circuit
// {
// private:
//     /* data */
//     int x;
// public:
//     Circuit(int){
//         this->x = x;
//     };
//     ~Circuit(){

//     };
//     void print(){
//         cout << "x = " << x << endl;
//     };
// };

Circuit::Circuit(int x)
{
    this->x = x;
}

Circuit::~Circuit()
{
}

void Circuit::print(){
    cout << "x = " << x << endl;
};

// int main(int argc, char const *argv[])
// {
//     /* code */
//     Circuit circuit = Circuit(1);
//     circuit.print();
//     return 0;
// }
