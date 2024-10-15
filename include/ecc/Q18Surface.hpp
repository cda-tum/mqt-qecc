#pragma once

#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <vector>
namespace ecc {
class Q18Surface : public Ecc {
public:
    Q18Surface(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Q18Surface,
               N_REDUNDANT_QUBITS,
               0,
               "Q18Surface",
               {{ANCILLA_WIDTH, "qeccX"}}},
              std::move(qc), measureFq) {}

    constexpr static QubitCount N_DATA_QUBITS    = 18;
    constexpr static QubitCount N_ANCILLA_QUBITS = 18;
    constexpr static QubitCount N_REDUNDANT_QUBITS =
            N_DATA_QUBITS + N_ANCILLA_QUBITS;
    constexpr static QubitCount ANCILLA_WIDTH = 8;

    constexpr static std::array<Qubit, N_DATA_QUBITS> DATA_QUBITS = {
            1, 3, 5, 6, 8, 10, 13, 15, 17, 18, 20, 22, 25, 27, 29, 30, 32, 34};
    constexpr static std::array<Qubit, N_ANCILLA_QUBITS> ANCILLA_INDICES = {
            0, 2, 4, 7, 9, 11, 12, 14, 16, 19, 21, 23, 24, 26, 28, 31, 33, 35};
    constexpr static std::array<Qubit, 4> ANCILLA_QUBITS_DECODE = {8, 13, 15, 20};
    constexpr static Qubit                X_INFORMATION         = 14;
    constexpr static std::array<Qubit, 3> LOGICAL_X             = {5, 10, 15};
    constexpr static std::array<Qubit, 3> LOGICAL_Z             = {20, 25, 30};

    //{a,{b,c}} == qubit a is checked by b and c
    std::map<std::size_t, std::vector<std::size_t>> qubitCorrectionX = {
            {1, {0, 2}}, {3, {2, 4}}, {5, {4}}, {6, {0, 12}}, {8, {2}}, {10, {4, 16}}, {13, {12}}, {15, {16}}, {18, {12, 24}}, {20, {26}}, {22, {16, 28}}, {25, {24, 26}}, {27, {26, 28}}, {29, {28}}, {30, {24}}};

    // temporarily deactivating phase correction
    /*std::map<std::size_t, std::vector<std::size_t>> qubitCorrectionZ = {
            {5, {11}},
            {6, {7}},
            {8, {7, 9}},
            {10, {9, 11}},
            {13, {7, 19}},
            {15, {9}},
            {17, {11, 23}},
            {20, {19}},
            {22, {23}},
            {25, {19, 31}},
            {27, {33}},
            {29, {23, 35}},
            {30, {31}},
            {32, {31, 33}},
            {34, {33, 35}}};*/

    static constexpr std::array<std::size_t, ANCILLA_WIDTH> X_CHECKS = {
            0, 2, 4, 12, 16, 24, 26, 28};
    // temporarily deactivating phase correction
    // static constexpr std::array<std::size_t, ANCILLA_WIDTH> zChecks = {7, 9,
    // 11, 19, 23, 31, 33, 35};

protected:
    void measureAndCorrect() override;

    void writeDecoding() override;

    void mapGate(const qc::Operation& gate) override;
};
} // namespace ecc
