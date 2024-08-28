#include "ecc/Q9Shor.hpp"

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
void Q9Shor::writeEncoding() {
    if (!isDecoded) {
        return;
    }
    isDecoded          = false;
    const auto nQubits = qcOriginal->getNqubits();
    for (Qubit i = 0; i < nQubits; i++) {
        std::array<qc::Control, 3> controls = {};
        for (std::size_t j = 0; j < controls.size(); j++) {
            controls.at(j) = {static_cast<Qubit>(i + 3 * j * nQubits),
                              qc::Control::Type::Pos};
            if (j > 0) {
                qcMapped->cx(controls[0], static_cast<Qubit>(i + 3 * j * nQubits));
            }
        }
        for (std::size_t j = 0; j < controls.size(); j++) {
            qcMapped->h(static_cast<Qubit>(i + 3 * j * nQubits));
            qcMapped->cx(controls.at(j),
                         static_cast<Qubit>(i + (3 * j + 1) * nQubits));
            qcMapped->cx(controls.at(j),
                         static_cast<Qubit>(i + (3 * j + 2) * nQubits));
        }
    }
    gatesWritten = true;
}

void Q9Shor::measureAndCorrect() {
    if (isDecoded) {
        return;
    }
    const auto nQubits = qcOriginal->getNqubits();
    const auto clStart = qcOriginal->getNcbits();
    for (Qubit i = 0; i < nQubits; i++) {
        // syntactic sugar for qubit indices
        std::array<Qubit, N_REDUNDANT_QUBITS>      qubits                  = {};
        std::array<Qubit, N_CORRECTING_BITS>       ancillaQubits           = {};
        std::array<qc::Control, N_CORRECTING_BITS> ancillaControls         = {};
        std::array<qc::Control, N_CORRECTING_BITS> negativeAncillaControls = {};
        for (std::size_t j = 0; j < qubits.size(); j++) {
            qubits.at(j) = static_cast<Qubit>(i + j * nQubits);
        }
        for (std::size_t j = 0; j < ancillaQubits.size(); j++) {
            ancillaQubits.at(j) =
                    static_cast<Qubit>(ecc.nRedundantQubits * nQubits + j);
            qcMapped->reset(ancillaQubits.at(j));
        }
        for (std::size_t j = 0; j < ancillaControls.size(); j++) {
            ancillaControls.at(j) = qc::Control{ancillaQubits.at(j)};
        }
        for (std::size_t j = 0; j < negativeAncillaControls.size(); j++) {
            negativeAncillaControls.at(j) =
                    qc::Control{ancillaQubits.at(j), qc::Control::Type::Neg};
        }

        // PREPARE measurements
        // --------------------------------------------------------
        for (Qubit const j : ancillaQubits) {
            qcMapped->h(j);
        }
        // x errors = indirectly via controlled z
        for (std::size_t j = 0; j < 3; j++) {
            qcMapped->cz(ancillaControls.at(2 * j), qubits.at(3 * j));
            qcMapped->cz(ancillaControls.at(2 * j), qubits.at(3 * j + 1));
            qcMapped->cz(ancillaControls.at(2 * j + 1), qubits.at(3 * j + 1));
            qcMapped->cz(ancillaControls.at(2 * j + 1), qubits.at(3 * j + 2));
        }

        // z errors = indirectly via controlled x/C-NOT
        for (std::size_t j = 0; j < 6; j++) {
            qcMapped->cx(ancillaControls[6], qubits.at(j));
            qcMapped->cx(ancillaControls[7], qubits.at(3 + j));
        }

        for (Qubit const j : ancillaQubits) {
            qcMapped->h(j);
        }

        // MEASURE ancilla qubits
        for (std::size_t j = 0; j < N_CORRECTING_BITS; j++) {
            qcMapped->measure(ancillaQubits.at(j), clStart + j);
        }

        // CORRECT
        // x, i.e. bit flip errors
        for (std::size_t j = 0; j < 3; j++) {
            const auto controlRegister =
                    std::make_pair(static_cast<Qubit>(clStart + 2 * j), 2);
            qcMapped->classicControlled(qc::X, qubits.at(3 * j), controlRegister, 1U);
            qcMapped->classicControlled(qc::X, qubits.at(3 * j + 2), controlRegister,
                                        2U);
            qcMapped->classicControlled(qc::X, qubits.at(3 * j + 1), controlRegister,
                                        3U);
        }

        // z, i.e. phase flip errors
        const auto controlRegister =
                std::make_pair(static_cast<Qubit>(clStart + 6), 2);
        qcMapped->classicControlled(qc::Z, qubits.at(0), controlRegister, 1U);
        qcMapped->classicControlled(qc::Z, qubits.at(6), controlRegister, 2U);
        qcMapped->classicControlled(qc::Z, qubits.at(3), controlRegister, 3U);
    }
}

void Q9Shor::writeDecoding() {
    if (isDecoded) {
        return;
    }
    const auto nQubits = qcOriginal->getNqubits();
    for (Qubit i = 0; i < nQubits; i++) {
        std::array<qc::Control, N_REDUNDANT_QUBITS> ci;
        for (Qubit j = 0; j < ci.size(); j++) {
            ci.at(j) = qc::Control{static_cast<Qubit>(i + j * nQubits),
                                   qc::Control::Type::Pos};
        }

        for (std::size_t j = 0; j < 3; j++) {
            std::array<Qubit, 3> targets = {
                    static_cast<Qubit>(i + 3 * j * nQubits),
                    static_cast<Qubit>(i + (3 * j + 1) * nQubits),
                    static_cast<Qubit>(i + (3 * j + 2) * nQubits)};
            qcMapped->cx(ci.at(3 * j), targets.at(1));
            qcMapped->cx(ci.at(3 * j), targets.at(2));
            qcMapped->mcx({ci.at(3 * j + 1), ci.at(3 * j + 2)}, targets.at(0));
            qcMapped->h(targets.at(0));
        }

        qcMapped->cx(ci[0], static_cast<Qubit>(i + 3 * nQubits));
        qcMapped->cx(ci[0], static_cast<Qubit>(i + 6 * nQubits));
        qcMapped->mcx({ci.at(3), ci.at(6)}, i);
    }
    isDecoded = true;
}

void Q9Shor::mapGate(const qc::Operation& gate) {
    if (isDecoded && gate.getType() != qc::Measure && gate.getType() != qc::H) {
        writeEncoding();
    }
    const auto nQubits = qcOriginal->getNqubits();
    auto       type    = qc::I;
    switch (gate.getType()) {
    case qc::I:
    case qc::Barrier:
        break;
    case qc::X:
        type = qc::Z;
        break;
    case qc::Y:
        type = qc::Y;
        break;
    case qc::Z:
        type = qc::X;
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
        return;
    default:
        gateNotAvailableError(gate);
    }
    for (std::size_t t = 0; t < gate.getNtargets(); t++) {
        auto i = gate.getTargets()[t];

        if (gate.getNcontrols() != 0U) {
            // Q9Shor code: put H gate before and after each control point, i.e. "cx
            // 0,1" becomes "h0; cz 0,1; h0"
            const auto& controls = gate.getControls();
            for (size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                qc::Controls controls2;
                for (const auto& ct : controls) {
                    controls2.insert(
                            qc::Control{static_cast<Qubit>(ct.qubit + j * nQubits), ct.type});
                    qcMapped->h(static_cast<Qubit>(ct.qubit + j * nQubits));
                }
                qcMapped->emplace_back<qc::StandardOperation>(
                        controls2, static_cast<Qubit>(i + j * nQubits), type);
                for (const auto& ct : controls) {
                    qcMapped->h(static_cast<Qubit>(ct.qubit + j * nQubits));
                }
            }
        } else {
            for (size_t j = 0; j < N_REDUNDANT_QUBITS; j++) {
                qcMapped->emplace_back<qc::StandardOperation>(
                        static_cast<Qubit>(i + j * nQubits), type);
            }
        }
    }
}
} // namespace ecc
