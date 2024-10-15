#pragma once

#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
namespace ecc {
class Q5Laflamme : public Ecc {
public:
    Q5Laflamme(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Q5Laflamme,
               N_REDUNDANT_QUBITS,
               N_CORRECTING_BITS,
               "Q5Laflamme",
               {{N_CORRECTING_BITS, "qecc"}, {1, "encode"}}},
              std::move(qc), measureFq) {}

protected:
    void writeEncoding() override;

    void measureAndCorrect() override;

    void writeDecoding() override;

    void mapGate(const qc::Operation& gate) override;

    static constexpr std::size_t N_REDUNDANT_QUBITS = 5;
    static constexpr std::size_t N_CORRECTING_BITS  = 4;

    static constexpr std::array<std::array<qc::OpType, N_REDUNDANT_QUBITS>, 4>
            STABILIZER_MATRIX = {{
                    {qc::X, qc::Z, qc::Z, qc::X, qc::I}, // c0
                    {qc::I, qc::X, qc::Z, qc::Z, qc::X}, // c1
                    {qc::X, qc::I, qc::X, qc::Z, qc::Z}, // c2
                    {qc::Z, qc::X, qc::I, qc::X, qc::Z}  // c3
            }};

    static constexpr std::array<Qubit, 8> DECODING_CORRECTION_VALUES = {
            1, 2, 4, 7, 8, 11, 13, 14}; // values with odd amount of '1' bits
};
} // namespace ecc
