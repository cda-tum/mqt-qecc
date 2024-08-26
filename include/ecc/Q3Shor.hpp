#pragma once

#include "Definitions.hpp"
#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <cstddef>
#include <memory>
#include <utility>
namespace ecc {
class Q3Shor : public Ecc {
public:
    Q3Shor(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Q3Shor, N_REDUNDANT_QUBITS, 2, "Q3Shor", {{2, "qecc"}}},
              std::move(qc), measureFq) {}

protected:
    void writeEncoding() override;

    void measureAndCorrect() override;

    void writeDecoding() override;

    void mapGate(const qc::Operation& gate) override;

    void addOperation(const qc::Controls& controls, const qc::Targets& targets,
                      qc::OpType type);

    static constexpr std::size_t N_REDUNDANT_QUBITS = 3;
};
} // namespace ecc
