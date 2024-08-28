#include "ecc/Q18Surface.hpp"

#include "ecc/Ecc.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/NonUnitaryOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
namespace ecc {
void Q18Surface::measureAndCorrect() {
    if (isDecoded) {
        return;
    }
    const auto nQubits    = qcOriginal->getNqubits();
    const auto clAncStart = qcOriginal->getNcbits();

    std::map<std::size_t, std::size_t> xCheckMasks;
    for (std::size_t j = 0; j < ANCILLA_WIDTH; j++) {
        xCheckMasks[X_CHECKS.at(j)] = static_cast<std::size_t>(1U) << j;
    }

    for (Qubit i = 0; i < nQubits; i++) {
        std::array<Qubit, N_REDUNDANT_QUBITS>       qubits        = {};
        std::array<qc::Control, N_REDUNDANT_QUBITS> controlQubits = {};
        for (std::size_t j = 0; j < qubits.size(); j++) {
            qubits.at(j) = static_cast<Qubit>(i + j * nQubits);
        }
        for (std::size_t j = 0; j < controlQubits.size(); j++) {
            controlQubits.at(j) = qc::Control{qubits.at(j)};
        }

        if (gatesWritten) {
            for (Qubit const ai : ANCILLA_INDICES) {
                qcMapped->reset(qubits.at(ai));
            }
        }

        // initialize ancillas: X-check
        for (const auto& [targetIndex, ancillaIndices] : qubitCorrectionX) {
            for (const auto ancilla : ancillaIndices) {
                qcMapped->cx(controlQubits.at(targetIndex), qubits.at(ancilla));
            }
        }

        // map ancillas to classical bit result
        for (std::size_t j = 0; j < X_CHECKS.size(); j++) {
            qcMapped->measure(qubits.at(X_CHECKS.at(j)), clAncStart + j);
        }

        // logic: classical control
        auto controlRegister =
                std::make_pair(static_cast<Qubit>(clAncStart), ANCILLA_WIDTH);
        for (const auto& [targetIndex, ancillaIndices] : qubitCorrectionX) {
            std::size_t mask = 0;
            for (std::size_t const ancillaIndex : ancillaIndices) {
                mask |= xCheckMasks[ancillaIndex];
            }
            qcMapped->classicControlled(qc::X, qubits.at(targetIndex),
                                        controlRegister, mask);
        }

        gatesWritten = true;
    }
}

void Q18Surface::writeDecoding() {
    if (isDecoded) {
        return;
    }
    const auto nQubits = qcOriginal->getNqubits();
    for (Qubit i = 0; i < nQubits; i++) {
        qcMapped->reset(static_cast<Qubit>(i + X_INFORMATION * nQubits));
        for (const Qubit qubit : ANCILLA_QUBITS_DECODE) {
            qcMapped->cx(static_cast<Qubit>(i + qubit * nQubits),
                         static_cast<Qubit>(i + X_INFORMATION * nQubits));
        }
        qcMapped->measure(static_cast<Qubit>(i + X_INFORMATION * nQubits), i);
        qcMapped->reset(i);
        qcMapped->cx(static_cast<Qubit>(i + X_INFORMATION * nQubits), i);
    }
    isDecoded = true;
}

void Q18Surface::mapGate(const qc::Operation& gate) {
    if (isDecoded && gate.getType() != qc::Measure) {
        writeEncoding();
    }
    const auto nQubits = qcOriginal->getNqubits();

    // no control gate decomposition is supported
    if (gate.isControlled() && gate.getType() != qc::Measure) {
        // multi-qubit gates are currently not supported
        gateNotAvailableError(gate);
    } else {
        static constexpr std::array<std::pair<Qubit, Qubit>, 6> SWAP_QUBIT_INDICES =
                {std::make_pair(1, 29), std::make_pair(3, 17), std::make_pair(6, 34),
                 std::make_pair(8, 22), std::make_pair(13, 27), std::make_pair(18, 32)};

        switch (gate.getType()) {
        case qc::Barrier:
        case qc::I:
            break;
        case qc::X:
            for (auto i : gate.getTargets()) {
                for (auto j : LOGICAL_X) {
                    qcMapped->x(static_cast<Qubit>(i + j * nQubits));
                }
            }
            break;
        case qc::H:
            // apply H gate to every data qubit
            // swap circuit along '/' axis
            for (auto i : gate.getTargets()) {
                for (const auto j : DATA_QUBITS) {
                    qcMapped->h(static_cast<Qubit>(i + j * nQubits));
                }
                for (auto pair : SWAP_QUBIT_INDICES) {
                    qcMapped->swap(
                            static_cast<Qubit>(i + static_cast<size_t>(pair.first) * nQubits),
                            static_cast<Qubit>(i +
                                               static_cast<size_t>(pair.second) * nQubits));
                }
                // qubits 5, 10, 15, 20, 25, 30 are along axis
            }
            break;
        case qc::Y:
            // Y = Z X
            for (auto i : gate.getTargets()) {
                for (auto j : LOGICAL_Z) {
                    qcMapped->z(static_cast<Qubit>(i + j * nQubits));
                }
                for (auto j : LOGICAL_X) {
                    qcMapped->x(static_cast<Qubit>(i + j * nQubits));
                }
            }
            break;
        case qc::Z:
            for (auto i : gate.getTargets()) {
                for (auto j : LOGICAL_Z) {
                    qcMapped->z(static_cast<Qubit>(i + j * nQubits));
                }
            }
            break;
        case qc::Measure:
            if (!isDecoded) {
                measureAndCorrect();
                writeDecoding();
            }
            if (const auto* measureGate =
                        dynamic_cast<const qc::NonUnitaryOperation*>(&gate)) {
                const auto& classics = measureGate->getClassics();
                const auto& targets  = measureGate->getTargets();
                for (std::size_t j = 0; j < classics.size(); j++) {
                    qcMapped->measure(targets.at(j), classics.at(j));
                }
            } else {
                throw std::runtime_error("Dynamic cast to NonUnitaryOperation failed.");
            }
            break;
        default:
            gateNotAvailableError(gate);
        }
    }
}
} // namespace ecc
