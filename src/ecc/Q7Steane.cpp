#include "ecc/Q7Steane.hpp"

#include "Definitions.hpp"
#include "ecc/Ecc.hpp"
#include "ir/operations/Control.hpp"
#include "ir/operations/NonUnitaryOperation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "ir/operations/StandardOperation.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
namespace ecc {
void Q7Steane::writeEncoding() {
    if (!isDecoded) {
        return;
    }
    isDecoded          = false;
    const auto nQubits = qcOriginal->getNqubits();
    // reset data qubits if necessary
    if (gatesWritten) {
        for (std::size_t i = 0; i < nQubits; i++) {
            for (std::size_t j = 1; j < N_REDUNDANT_QUBITS; j++) {
                qcMapped->reset(static_cast<Qubit>(i + j * nQubits));
            }
        }
    }
    measureAndCorrectSingle(true);
}

void Q7Steane::measureAndCorrect() {
    if (isDecoded) {
        return;
    }
    measureAndCorrectSingle(true);
    measureAndCorrectSingle(false);
}

void Q7Steane::measureAndCorrectSingle(bool xSyndrome) {
    const auto nQubits    = qcOriginal->getNqubits();
    const auto ancStart   = nQubits * ecc.nRedundantQubits;
    const auto clAncStart = qcOriginal->getNcbits();
    const auto controlRegister =
            std::make_pair(static_cast<Qubit>(clAncStart), N_CORRECTING_BITS);

    for (Qubit i = 0; i < nQubits; i++) {
        if (gatesWritten) {
            for (std::size_t j = 0; j < ecc.nCorrectingBits; j++) {
                qcMapped->reset(static_cast<Qubit>(ancStart + j));
            }
        }

        std::array<qc::Control, 3> controls = {};
        for (std::size_t j = 0; j < ecc.nCorrectingBits; j++) {
            qcMapped->h(static_cast<Qubit>(ancStart + j));
            controls.at(j) = qc::Control{static_cast<Qubit>(ancStart + j)};
        }

        // K1: UIUIUIU
        // K2: IUUIIUU
        // K3: IIIUUUU
        for (std::size_t c = 0; c < controls.size(); c++) {
            for (std::size_t q = 0; q < ecc.nRedundantQubits; q++) {
                if (((q + 1) & (static_cast<std::size_t>(1U) << c)) != 0) {
                    const auto target = static_cast<Qubit>(i + nQubits * q);
                    if (xSyndrome) {
                        qcMapped->cx(controls.at(c), target);
                    } else {
                        qcMapped->cz(controls.at(c), target);
                    }
                }
            }
        }

        for (std::size_t j = 0; j < ecc.nCorrectingBits; j++) {
            qcMapped->h(static_cast<Qubit>(ancStart + j));
            qcMapped->measure(static_cast<Qubit>(ancStart + j), clAncStart + j);
        }

        // correct Z_i for i+1 = c0*1+c1*2+c2*4
        // correct X_i for i+1 = c3*1+c4*2+c5*4
        const auto opType = xSyndrome ? qc::Z : qc::X;
        for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
            qcMapped->classicControlled(opType, static_cast<Qubit>(i + j * nQubits),
                                        controlRegister, j + 1U);
        }
        gatesWritten = true;
    }
}

void Q7Steane::writeDecoding() {
    if (isDecoded) {
        return;
    }
    const auto nQubits    = qcOriginal->getNqubits();
    const auto clAncStart = qcOriginal->getNcbits();

    // use exiting registers qeccX and qeccZ for decoding
    const auto controlRegister =
            std::make_pair(static_cast<Qubit>(clAncStart), N_CORRECTING_BITS);

    for (Qubit i = 0; i < nQubits; i++) {
        // #|###|###
        // 0|111|111
        // odd amount of 1's -> x.at(0) = 1
        // measure from index 1 (not 0) to 6, =qubit 2 to 7

        qcMapped->measure(static_cast<Qubit>(i + 1 * nQubits), clAncStart);
        qcMapped->measure(static_cast<Qubit>(i + 2 * nQubits), clAncStart + 1);
        qcMapped->measure(static_cast<Qubit>(i + 3 * nQubits), clAncStart + 2);
        for (auto value : DECODING_CORRECTION_VALUES) {
            qcMapped->classicControlled(qc::X, i, controlRegister, value);
        }
        qcMapped->measure(static_cast<Qubit>(i + 4 * nQubits), clAncStart);
        qcMapped->measure(static_cast<Qubit>(i + 5 * nQubits), clAncStart + 1);
        qcMapped->measure(static_cast<Qubit>(i + 6 * nQubits), clAncStart + 2);
        for (auto value : DECODING_CORRECTION_VALUES) {
            qcMapped->classicControlled(qc::X, i, controlRegister, value);
        }
    }
    isDecoded = true;
}

void Q7Steane::addSOperation(const qc::Controls& controls,
                             const qc::Targets&  targets,
                             const qc::OpType    type) {
    const auto numTargets  = targets.size();
    const auto numControls = controls.size();
    const auto numQubits   = qcOriginal->getNqubits();

    for (std::size_t k = 0; k < numTargets; k++) {
        auto i = targets[k];
        if (numControls > 0) {
            for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                qc::Controls controls2;
                for (const auto& ct : controls) {
                    controls2.emplace(static_cast<Qubit>(ct.qubit + j * numQubits),
                                      ct.type);
                }
                qcMapped->emplace_back<qc::StandardOperation>(
                        controls2, static_cast<Qubit>(i + j * numQubits), type);
                qcMapped->emplace_back<qc::StandardOperation>(
                        controls2, static_cast<Qubit>(i + j * numQubits), type);
                qcMapped->emplace_back<qc::StandardOperation>(
                        controls2, static_cast<Qubit>(i + j * numQubits), type);
            }
        } else {
            for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                qcMapped->emplace_back<qc::StandardOperation>(
                        static_cast<Qubit>(i + j * numQubits), type);
                qcMapped->emplace_back<qc::StandardOperation>(
                        static_cast<Qubit>(i + j * numQubits), type);
                qcMapped->emplace_back<qc::StandardOperation>(
                        static_cast<Qubit>(i + j * numQubits), type);
            }
        }
    }
}

void Q7Steane::mapGate(const qc::Operation& gate) {
    if (isDecoded && gate.getType() != qc::Measure) {
        writeEncoding();
    }
    const QubitCount nQubits = qcOriginal->getNqubits();
    switch (gate.getType()) {
    case qc::I:
    case qc::Barrier:
        break;
    case qc::X:
    case qc::H:
    case qc::Y:
    case qc::Z:
        for (auto i : gate.getTargets()) {
            if (gate.getNcontrols() != 0U) {
                const auto& ctrls = gate.getControls();
                for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                    qc::Controls ctrls2;
                    for (const auto& ct : ctrls) {
                        ctrls2.insert(qc::Control{
                                static_cast<Qubit>(ct.qubit + j * nQubits), ct.type});
                    }
                    qcMapped->emplace_back<qc::StandardOperation>(
                            ctrls2, static_cast<Qubit>(i + j * nQubits), gate.getType());
                }
            } else {
                for (std::size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                    qcMapped->emplace_back<qc::StandardOperation>(
                            static_cast<Qubit>(i + j * nQubits), gate.getType());
                }
            }
        }
        break;
        // locigal S = 3 physical S's
    case qc::S:
    case qc::Sdg:
        addSOperation(gate.getControls(), gate.getTargets(), gate.getType());
        break;
    case qc::Measure:
        if (!isDecoded) {
            measureAndCorrect();
            writeDecoding();
        }
        if (const auto* measureGate =
                    dynamic_cast<const qc::NonUnitaryOperation*>(&gate)) {
            for (std::size_t j = 0; j < measureGate->getNclassics(); j++) {
                auto classicalRegisterName =
                        qcOriginal->getClassicalRegister(measureGate->getTargets().at(j));
                if (!classicalRegisterName.empty()) {
                    qcMapped->measure(
                            static_cast<Qubit>(measureGate->getClassics().at(j)),
                            {classicalRegisterName, measureGate->getTargets().at(j)});
                } else {
                    qcMapped->measure(
                            static_cast<Qubit>(measureGate->getClassics().at(j)),
                            measureGate->getTargets().at(j));
                }
            }
        } else {
            throw std::runtime_error("Dynamic cast to NonUnitaryOperation failed.");
        }

        break;
    case qc::T:
    case qc::Tdg:
        for (auto i : gate.getTargets()) {
            if (gate.getControls().empty()) {
                qcMapped->cx(static_cast<Qubit>(i + 6 * nQubits),
                             static_cast<Qubit>(i + 5 * nQubits));
                qcMapped->cx(static_cast<Qubit>(i + 5 * nQubits),
                             static_cast<Qubit>(i + 0 * nQubits));
                if (gate.getType() == qc::T) {
                    qcMapped->t(static_cast<Qubit>(i + 0 * nQubits));
                } else {
                    qcMapped->tdg(static_cast<Qubit>(i + 0 * nQubits));
                }
                qcMapped->cx(static_cast<Qubit>(i + 5 * nQubits),
                             static_cast<Qubit>(i + 0 * nQubits));
                qcMapped->cx(static_cast<Qubit>(i + 6 * nQubits),
                             static_cast<Qubit>(i + 5 * nQubits));
            } else {
                gateNotAvailableError(gate);
            }
        }
        break;
    default:
        gateNotAvailableError(gate);
    }
}
} // namespace ecc
