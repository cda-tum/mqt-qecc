#pragma once

#include "Definitions.hpp"
#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
namespace ecc {
class Q7Steane : public Ecc {
public:
    Q7Steane(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Q7Steane,
               N_REDUNDANT_QUBITS,
               N_CORRECTING_BITS,
               "Q7Steane",
               {{N_CORRECTING_BITS, "qecc"}}},
              std::move(qc), measureFq) {}

protected:
    void writeEncoding() override;

    void measureAndCorrect() override;
    void measureAndCorrectSingle(bool xSyndrome);

    void writeDecoding() override;

    void mapGate(const qc::Operation& gate) override;

    void addSOperation(const qc::Controls& controls, const qc::Targets& targets,
                       qc::OpType type);

    static constexpr std::size_t N_REDUNDANT_QUBITS = 7;
    static constexpr std::size_t N_CORRECTING_BITS  = 3;

    static constexpr std::array<Qubit, 4> DECODING_CORRECTION_VALUES = {
            1, 2, 4, 7}; // values with odd amount of '1' bits
};
} // namespace ecc
