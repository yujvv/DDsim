/*
 * This file is part of the MQT DD Package which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/quantum_dd/ for more information.
 */

#include "dd/Export.hpp"
#include "dd/GateMatrixDefinitions.hpp"
#include "dd/Package.hpp"
#include "dd/Edge.hpp"

#include <complex>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <cmath>
#include <stack>
#include <omp.h>

using namespace dd::literals;

// auto bellCicuit1(std::unique_ptr<dd::Package<>>& dd) {
//     /***** define Hadamard gate acting on q0 *****/
//     auto hGate = dd->makeGateDD(dd::Hmat, 2, 0);

//     /***** define cx gate with control q0 and target q1 *****/
//     auto cxGate = dd->makeGateDD(dd::Xmat, 2, 0_pc, 1);

//     //Multiply matrices to get functionality of circuit
//     return dd->multiply(cxGate, hGate);
// }

// auto bellCicuit2(std::unique_ptr<dd::Package<>>& dd) {
//     /***** define Hadamard gate acting on q1 *****/
//     auto hGateQ1 = dd->makeGateDD(dd::Hmat, 2, 1);

//     /***** define Hadamard gate acting on q0 *****/
//     auto hGateQ0 = dd->makeGateDD(dd::Hmat, 2, 0);

//     /***** define cx gate with control q1 and target q0 *****/
//     auto cxGate = dd->makeGateDD(dd::Xmat, 2, 1_pc, 0);

//     //Multiply matrices to get functionality of circuit
//     return dd->multiply(dd->multiply(hGateQ1, hGateQ0), dd->multiply(cxGate, hGateQ1));
// }

// int main() {                           // NOLINT(bugprone-exception-escape)
//     dd::Package<>::printInformation(); // uncomment to print various sizes of structs and arrays
//     //Initialize package
//     auto dd = std::make_unique<dd::Package<>>(4);

//     // create Bell circuit 1
//     auto bellCircuit1 = bellCicuit1(dd);

//     // create Bell circuit 2
//     auto bellCircuit2 = bellCicuit2(dd);

//     /***** Equivalence checking *****/
//     if (bellCircuit1 == bellCircuit2) {
//         std::cout << "Circuits are equal!" << std::endl;
//         std::cout << "debug" << std::endl;
//     } else {
//         std::cout << "Circuits are not equal!" << std::endl;
//     }

//     /***** Simulation *****/
//     //Generate vector in basis state |00>
//     auto zeroState = dd->makeZeroState(2);

//     //Simulate the bell_circuit with initial state |00>
//     auto bellState  = dd->multiply(bellCircuit1, zeroState);
//     auto bellState2 = dd->multiply(bellCircuit2, zeroState);

//     //print result
//     dd->printVector(bellState);

//     std::cout << "Bell states have a fidelity of " << dd->fidelity(bellState, bellState2) << "\n";
//     std::cout << "Bell state and zero state have a fidelity of " << dd->fidelity(bellState, zeroState) << "\n";

//     /***** Custom gates *****/
//     // define, e.g., Pauli-Z matrix
//     dd::GateMatrix m;
//     m[0] = {1., 0.};
//     m[1] = {0., 0.};
//     m[2] = {0., 0.};
//     m[3] = {-1., 0.};

//     auto myZGate = dd->makeGateDD(m, 1, 0);
//     std::cout << "DD of my gate has size " << dd->size(myZGate) << std::endl;

//     // compute (partial) traces
//     auto partTrace = dd->partialTrace(dd->makeIdent(2), {true, true});
//     auto fullTrace = dd->trace(dd->makeIdent(4));
//     std::cout << "Identity function for 4 qubits has trace: " << fullTrace << std::endl;

//     /***** print DDs as SVG file *****/
//     dd::export2Dot(bellCircuit1, "bell_circuit1.dot", false);
//     dd::export2Dot(bellCircuit2, "bell_circuit2.dot");
//     dd::export2Dot(bellState, "bell_state.dot", true);
//     dd::export2Dot(partTrace, "partial_trace.dot");

//     /***** print statistics *****/
//     dd->statistics();
//     dd->garbageCollect(true);
// }

double randDouble(double min, double max) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(min, max);
  return dis(gen);
}

std::complex<double> randComplex() {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-1.0, 1.0);
  return std::complex<double>(dis(gen), dis(gen));
}

auto randomState(std::unique_ptr<dd::Package<>>& dd, int num_qubits) {

    std::vector<dd::BasisStates> state(num_qubits);
    std::mt19937 gen(time(nullptr));
    std::uniform_int_distribution<> dis(0, 1);

    // auto basis_state = dd->makeBasisState(static_cast<dd::QubitCount>(num_qubits));
    for (int i = 0; i < num_qubits; i++) {
        state[i] = dis(gen) ? dd::BasisStates::one : dd::BasisStates::zero;
    }

    dd::QubitCount num = num_qubits;

    
    auto basis_state = dd->makeBasisState(num, state);

    return basis_state;
}


auto multiply_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    std::cout << "The Test of multiplying hGate by cxGate ---" << std::endl;

    // std::cout << "Input The dimension of matrix:";
    // int d_matrix;
    // std::cin >> d_matrix;

    /***** define Hadamard gate acting on q0 *****/
    auto hGate = dd->makeGateDD(dd::Hmat, num_qubits, 0);

    /***** define cx gate with control q0 and target q1 *****/
    auto cxGate = dd->makeGateDD(dd::Xmat, num_qubits, 0_pc, 1);

    auto res = dd->multiply(cxGate, hGate);
    
    // std::cout << "Please enter the Number of calculations:";
    // int count;
    // std::cin >> count;

    auto t1 = std::chrono::high_resolution_clock::now();  
    for (int i = 0; i < num_try; i++) {
        res = dd->multiply(cxGate, res);
    }
    //Multiply matrices to get functionality of circuit
    auto t2 = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double, std::milli> ms = t2 - t1;

    std::cout << "time is  " <<ms.count()<<" ms"<< std::endl;
    return ms.count();
}

auto RD_multiply_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    using namespace dd;
        std::cout << "The Test of random multiplying Gate by Gate ---" << std::endl;

        GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

        std::vector<mEdge> gate_queue;
        const std::size_t NGATES = num_try;
        const uint64_t NQUBITS = num_qubits;

        std::mt19937_64 rng;
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(ss);

        std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
        std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );


        auto t1 = std::chrono::high_resolution_clock::now();
        for(std::size_t i = 0; i < NGATES; i++){
            gate_queue.push_back(dd->makeGateDD(gates[gate_dist(rng)], NQUBITS, qubit_dist(rng)));    
        }
        auto result = dd->multiply(gate_queue[0], gate_queue[1]);
        for(auto i = 2; i < NGATES; i++) result = dd->multiply(result, gate_queue[i]);
        auto t2 = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout<< "time is " <<ms.count()<<" ms"<<std::endl;

        return ms.count();
}

auto RD_add_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    using namespace dd;
        std::cout << "The Test of random addtion Gate by Gate ---" << std::endl;

        GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

        std::vector<mEdge> gate_queue;
        const std::size_t NGATES = num_try;
        const uint64_t NQUBITS = num_qubits;

        std::mt19937_64 rng;
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(ss);

        std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
        std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );


        auto t1 = std::chrono::high_resolution_clock::now();
        for(std::size_t i = 0; i < NGATES; i++){
            gate_queue.push_back(dd->makeGateDD(gates[gate_dist(rng)], NQUBITS, qubit_dist(rng)));    
        }
        auto result = dd->add(gate_queue[0], gate_queue[1]);
        for(auto i = 2; i < NGATES; i++) result = dd->add(result, gate_queue[i]);
        auto t2 = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout<< "time is " <<ms.count()<<" ms"<<std::endl;

        return ms.count();
}

auto vector_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    std::cout << "The Test of multiplying hGate by vector ---" << std::endl;

    // std::cout << "Input The dimension of vector and matrix:";
    // int d_matrix;
    // std::cin >> d_matrix;

    // Generate vector in basis state |00>
    auto zeroState = dd->makeZeroState(num_qubits);
    auto hGate = dd->makeGateDD(dd::Hmat, num_qubits, 0);

    auto res = dd->multiply(hGate, zeroState);
    
    // std::cout << "Please enter the Number of calculations:";
    // int count;
    // std::cin >> count;

    auto t1 = std::chrono::high_resolution_clock::now(); 
    for (int i = 0; i < num_try; i++) {
        res = dd->multiply(hGate, res);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double, std::milli> ms = t2 - t1;

    std::cout << "time is  " <<ms.count()<<" ms"<< std::endl;

    // end = clock();

    // std::cout << "time is  " << std::fixed << std::setprecision(10) << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
    return ms.count();
}

auto PD_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    std::cout << "The Test of Probability distribution of quantum states ---" << std::endl;
    auto RDState = randomState(dd, num_qubits);

    auto t1 = std::chrono::high_resolution_clock::now(); 
    for (int i = 0; i < num_try; i++) {
        auto probabilities = dd->determineMeasurementProbabilities(RDState, 0, true);
        // myStack.pop();
        // double p0 = probabilities[0];
    }

    auto t2 = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double, std::milli> ms = t2 - t1;

    std::cout << "time is  " <<ms.count()<<" ms"<< std::endl;
    
    return ms.count();
}

auto innerProduct_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    std::cout << "The Test of Inner Product of quantum states ---" << std::endl;
    clock_t start,end;

    auto RDState1 = randomState(dd, num_qubits);
    auto RDState2 = randomState(dd, num_qubits);
    // auto RDState = dd->makeZeroState(2);

    auto t1 = std::chrono::high_resolution_clock::now();  
    for (int i = 0; i < num_try; i++) {
        auto probabilities = dd->innerProduct(RDState1, RDState2);
        // myStack.pop();
        // double p0 = probabilities[0];
    }

    auto t2 = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double, std::milli> ms = t2 - t1;

    std::cout << "time is  " <<ms.count()<<" ms"<< std::endl;
    
    return ms.count();
}

auto RD_multiply_mv_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
    using namespace dd;
        std::cout << "The Test of multiplying Gate by state ---" << std::endl;
        
        auto RDState = randomState(dd, num_qubits);

        GateMatrix gates[] = {Imat, Hmat, Xmat, Ymat, Zmat, Smat, Sdagmat, Tmat, Tdagmat, SXmat, SXdagmat, Vmat, Vdagmat};

        std::vector<mEdge> gate_queue;
        const std::size_t NGATES = num_try;
        const uint64_t NQUBITS = num_qubits;

        std::mt19937_64 rng;
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(ss);

        std::uniform_int_distribution<int> qubit_dist(0, NQUBITS-1);
        std::uniform_int_distribution<int> gate_dist(0, gates->size() - 1 );

        auto t1 = std::chrono::high_resolution_clock::now();
        for(std::size_t i = 0; i < NGATES; i++){
            gate_queue.push_back(dd->makeGateDD(gates[gate_dist(rng)], NQUBITS, qubit_dist(rng)));    
        }
        auto result = dd->multiply(gate_queue[0], RDState);
        for(auto i = 2; i < NGATES; i++) result = dd->multiply(gate_queue[i], result);
        auto t2 = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout<< "time is " <<ms.count()<<" ms"<<std::endl;

        return ms.count();
}


// auto normalize_test(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
//     std::cout << "The Test of Normalize of Vector ---" << std::endl;
//     clock_t start,end;

//     std::stack<Edge<Node>> myStack;

//     for (int i = 0; i < num_try; i++) {
//         Edge<Node> e = {0, num_qubits - 1, randDouble(0.0, 1.0)};
//         myStack.push(e);
//     }

//     start = clock();  
//     for (int i = 0; i < num_try; i++) {
//         Edge<Node> e = myStack.top();
//         Edge<Node> normalized = e.normalize(e, false);
//         myStack.pop();
//         // double p0 = probabilities[0];
//     }

//     end = clock();

//     std::cout << "time is  " << std::fixed << std::setprecision(10) << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
//     return end-start;
// }


// auto randomGateMatrix(std::unique_ptr<dd::Package<>>& dd, int num_qubits) {
//     dd::vEdge zeroState = dd.makeZeroState(num_qubits);
//     std::vector<std::vector<std::complex<double>>> mat(3, std::vector<std::complex<double>>(3));
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             mat[i][j] = randComplex();
//         }
//     }
//     // Convert the matrix to a DD object  todo
//     auto gate = dd->makeDiagonal(mat);
//     // Apply the quantum gate to the zero state
//     dd::vEdge result = dd.apply(gate, zeroState);

//     return gate;
// }



// auto random_Matrix_multiply(std::unique_ptr<dd::Package<>>& dd, int num_qubits, int num_try) {
//     std::cout << "The Test of multiplying random Gate by random Gate ---" << std::endl;
//     clock_t start,end;

//     auto rndMatrix1 = randomGateMatrix(dd, num_qubits);
//     auto rndMatrix2 = randomGateMatrix(dd, num_qubits);
//     // auto rndState = randomState(dd, num_qubits);
//     auto res = dd->multiply(rndMatrix1, rndMatrix2);

//     start = clock();
//     for (int i = 0; i < num_try; i++) {
//         res = dd->multiply(rndMatrix1, res);
//     }
//     end = clock();

//     std::cout << "time is  " << std::fixed << std::setprecision(10) << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
//     return end-start;
// }



int main() { 
    auto dd = std::make_unique<dd::Package<>>(20);

    int num_qubits, num_try;
    std::cout << "Input the Number of qubits: ";
    std::cin >> num_qubits;
    std::cout << "Input the Number of attempts: ";
    std::cin >> num_try;
 

    // // Multiplication of quantum gates
    // auto time_multiply_mm = multiply_test(dd, num_qubits, num_try);
    // Multiplication of Quantum Gate and Quantum State
    auto time_multiply_mv = vector_test(dd, num_qubits, num_try);
    // Calculation of probability distribution of quantum states
    auto time_PD = PD_test(dd, num_qubits, num_try);
    // Calculation of Inner Product (Fidelity) of quantum states
    auto time_innerProduct = innerProduct_test(dd, num_qubits, num_try);
    // Multiplication of random quantum gates
    auto time_RD_multiply_mm = RD_multiply_test(dd, num_qubits, num_try);
    // Addtion of random quantum gates
    auto time_RD_add_mm = RD_add_test(dd, num_qubits, num_try);
    // Multiplication of random quantum Gates by quantum State
    auto time_RD_multiply_mv = RD_multiply_mv_test(dd, num_qubits, num_try);


    // Calculation of Normalize of Vector
    // auto time_normalize = normalize_test(dd, num_qubits, num_try);


    // Edge add(const Edge& x, const Edge& y) 

    // 计算4维度的identity的迹
    // auto fullTrace = dd->trace(dd->makeIdent(4));
    
    //Matrix (conjugate) transpose

    // test
    // #pragma omp parallel for num_threads(4)
    // for (int i = 0; i < 10; i++) {
    //     std::cout << "Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
    // }

}