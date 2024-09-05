#pragma once

#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <set>
#include <utility>

// Reference to this ecc in https://arxiv.org/pdf/1608.05053.pdf
namespace ecc {
class Q9Surface : public Ecc {
public:
    Q9Surface(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Q9Surface,
               N_REDUNDANT_QUBITS,
               N_CORRECTING_BITS,
               "Q9Surface",
               {{ANCILLA_WIDTH, "qeccX"}, {ANCILLA_WIDTH, "qeccZ"}}},
              std::move(qc), measureFq) {}

protected:
    void measureAndCorrect() override;

    void writeDecoding() override;

    void mapGate(const qc::Operation& gate) override;

private:
    //{a,{b,c}} == qubit a is checked by b and c
    static constexpr std::size_t N_REDUNDANT_QUBITS = 9;
    static constexpr std::size_t N_CORRECTING_BITS  = 8;
    static constexpr std::size_t ANCILLA_WIDTH      = 4;

    std::array<std::set<std::size_t>, N_REDUNDANT_QUBITS> qubitCorrectionX = {
            {{1}, {3}, {3}, {1, 4}, {3, 4}, {3, 6}, {4}, {4}, {6}}};
    std::array<std::set<std::size_t>, N_REDUNDANT_QUBITS> qubitCorrectionZ = {
            {{2}, {0, 2}, {0}, {2}, {2, 5}, {5}, {7}, {5, 7}, {5}}};
    static constexpr std::array<Qubit, ANCILLA_WIDTH>       X_ANCILLA_QUBITS   = {1, 3, 4,
                                                                                  6};
    static constexpr std::array<Qubit, ANCILLA_WIDTH>       Z_ANCILLA_QUBITS   = {0, 2, 5,
                                                                                  7};
    std::set<Qubit>                                         uncorrectedXQubits = {2, 6};
    std::set<Qubit>                                         uncorrectedZQubits = {0, 8};
    static constexpr std::array<Qubit, 3>                   LOGICAL_X          = {2, 4, 6};
    static constexpr std::array<Qubit, 3>                   LOGICAL_Z          = {0, 4, 8};
    static constexpr std::array<std::pair<Qubit, Qubit>, 4> SWAP_INDICES       = {
            {{0, 6}, {3, 7}, {2, 8}, {1, 5}}};
};
} // namespace ecc
