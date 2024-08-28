#include "ecc/Q9Surface.hpp"

#include "ecc/Ecc.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/NonUnitaryOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
namespace ecc {
void Q9Surface::measureAndCorrect() {
    if (isDecoded) {
        return;
    }
    const auto nQubits    = qcOriginal->getNqubits();
    const auto ancStart   = qcOriginal->getNqubits() * ecc.nRedundantQubits;
    const auto clAncStart = qcOriginal->getNcbits();
    for (std::size_t i = 0; i < nQubits; i++) {
        std::array<Qubit, N_REDUNDANT_QUBITS>       qubits          = {};
        std::array<Qubit, N_CORRECTING_BITS>        ancillaQubits   = {};
        std::array<qc::Control, N_CORRECTING_BITS>  ancillaControls = {};
        std::array<qc::Control, N_REDUNDANT_QUBITS> controlQubits   = {};
        for (std::size_t j = 0; j < qubits.size(); j++) {
            qubits.at(j) = static_cast<Qubit>(i + j * nQubits);
        }
        for (std::size_t j = 0; j < ancillaQubits.size(); j++) {
            ancillaQubits.at(j) = static_cast<Qubit>(ancStart + j);
        }
        if (gatesWritten) {
            for (const auto ancillaQubit : ancillaQubits) {
                qcMapped->reset(ancillaQubit);
            }
        }
        for (std::size_t j = 0; j < ancillaControls.size(); j++) {
            ancillaControls.at(j) = qc::Control{ancillaQubits.at(j)};
        }
        for (std::size_t j = 0; j < controlQubits.size(); j++) {
            controlQubits.at(j) = qc::Control{qubits.at(j)};
        }

        // X-type check (z error) on a0, a2, a5, a7: cx ancillaQubits->qubits
        // Z-type check (x error) on a1, a3, a4, a6: cz ancillaQubits->qubits = cx
        // qubits->ancillaQubits, no hadamard gate
        for (auto q : Z_ANCILLA_QUBITS) {
            qcMapped->h(ancillaQubits.at(q));
        }

        for (std::size_t q = 0; q < qubitCorrectionZ.size(); q++) {
            for (auto c : qubitCorrectionZ.at(q)) {
                qcMapped->cx(ancillaControls.at(c), qubits.at(q));
            }
        }

        for (std::size_t q = 0; q < qubitCorrectionX.size(); q++) {
            for (auto c : qubitCorrectionX.at(q)) {
                qcMapped->cx(controlQubits.at(q), ancillaQubits.at(c));
            }
        }

        for (std::size_t j = 0; j < Z_ANCILLA_QUBITS.size(); j++) {
            qcMapped->h(ancillaQubits.at(Z_ANCILLA_QUBITS.at(j)));
            qcMapped->measure(ancillaQubits.at(Z_ANCILLA_QUBITS.at(j)),
                              clAncStart + j);
            qcMapped->measure(ancillaQubits.at(X_ANCILLA_QUBITS.at(j)),
                              clAncStart + 4 + j);
        }

        // correction
        auto controlRegister =
                std::make_pair(static_cast<Qubit>(clAncStart), ANCILLA_WIDTH);
        for (std::size_t q = 0; q < qubitCorrectionZ.size(); q++) {
            if (uncorrectedZQubits.count(static_cast<Qubit>(q)) == 0) {
                std::size_t mask = 0;
                for (std::size_t c = 0; c < Z_ANCILLA_QUBITS.size(); c++) {
                    if (qubitCorrectionZ.at(q).count(Z_ANCILLA_QUBITS.at(c)) > 0) {
                        mask |= (static_cast<std::size_t>(1U) << c);
                    }
                }
                qcMapped->classicControlled(qc::Z, qubits.at(q), controlRegister, mask);
            }
        }
        controlRegister = std::make_pair(
                static_cast<Qubit>(clAncStart + ANCILLA_WIDTH), ANCILLA_WIDTH);
        for (std::size_t q = 0; q < qubitCorrectionX.size(); q++) {
            if (uncorrectedXQubits.count(static_cast<Qubit>(q)) == 0) {
                std::size_t mask = 0;
                for (std::size_t c = 0; c < X_ANCILLA_QUBITS.size(); c++) {
                    if (qubitCorrectionX.at(q).count(X_ANCILLA_QUBITS.at(c)) > 0) {
                        mask |= (static_cast<std::size_t>(1U) << c);
                    }
                }
                qcMapped->classicControlled(qc::X, qubits.at(q), controlRegister, mask);
            }
        }

        gatesWritten = true;
    }
}

void Q9Surface::writeDecoding() {
    if (isDecoded) {
        return;
    }
    const auto nQubits = qcOriginal->getNqubits();
    for (std::size_t i = 0; i < nQubits; i++) {
        // measure 0, 4, 8. state = m0*m4*m8
        qcMapped->measure(static_cast<Qubit>(i), i);
        qcMapped->measure(static_cast<Qubit>(i + 4 * nQubits), i);
        qcMapped->measure(static_cast<Qubit>(i + 8 * nQubits), i);
        qcMapped->cx(static_cast<Qubit>(i + 4 * nQubits), static_cast<Qubit>(i));
        qcMapped->cx(static_cast<Qubit>(i + 8 * nQubits), static_cast<Qubit>(i));
        qcMapped->measure(static_cast<Qubit>(i), i);
    }
    isDecoded = true;
}

void Q9Surface::mapGate(const qc::Operation& gate) {
    if (isDecoded && gate.getType() != qc::Measure) {
        writeEncoding();
    }
    const auto nQubits = qcOriginal->getNqubits();

    if ((gate.getNcontrols() != 0U) && gate.getType() != qc::Measure) {
        gateNotAvailableError(gate);
    }

    switch (gate.getType()) {
    case qc::I:
    case qc::Barrier:
        break;
    case qc::X:
        for (auto i : gate.getTargets()) {
            for (auto j : LOGICAL_X) {
                qcMapped->x(static_cast<Qubit>(i + j * nQubits));
            }
        }
        break;
    case qc::H:
        for (auto i : gate.getTargets()) {
            for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                qcMapped->h(static_cast<Qubit>(i + j * nQubits));
            }
            for (auto pair : SWAP_INDICES) {
                qcMapped->swap(static_cast<Qubit>(i + pair.first * nQubits),
                               static_cast<Qubit>(i + pair.second * nQubits));
            }
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
            for (std::size_t j = 0; j < measureGate->getNclassics(); j++) {
                qcMapped->measure(measureGate->getTargets().at(j),
                                  measureGate->getClassics().at(j));
            }
        } else {
            throw std::runtime_error("Dynamic cast to NonUnitaryOperation failed.");
        }
        break;
    default:
        gateNotAvailableError(gate);
    }
}
} // namespace ecc
