#pragma once

#include "Ecc.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"

#include <cstddef>
#include <memory>
#include <utility>
namespace ecc {
class Id : public Ecc {
public:
    Id(std::shared_ptr<qc::QuantumComputation> qc, std::size_t measureFq)
        : Ecc({ID::Id, 1, 0, "Id", {}}, std::move(qc), measureFq) {}

protected:
    void writeEncoding() override {};

    void measureAndCorrect() override {};

    void writeDecoding() override {};

    void mapGate(const qc::Operation& gate) override {
        qcMapped->emplace_back(gate.clone());
    }
};
} // namespace ecc
