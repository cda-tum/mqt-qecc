#include "ecc/Ecc.hpp"

#include "ir/QuantumComputation.hpp"

#include <cstddef>
#include <memory>

namespace ecc {
void Ecc::initMappedCircuit() {
    qcOriginal->stripIdleQubits(true, false);
    qcMapped->addQubitRegister(getNOutputQubits(qcOriginal->getNqubits()));
    auto cRegs = qcOriginal->getCregs();
    for (auto const& [regName, regBits] : cRegs) {
        qcMapped->addClassicalRegister(regBits.second, regName);
    }
    for (auto const& [regBits, regName] : ecc.classicalRegisters) {
        qcMapped->addClassicalRegister(regBits, regName);
    }
}

std::shared_ptr<qc::QuantumComputation> Ecc::apply() {
    initMappedCircuit();

    writeEncoding();
    isDecoded = false;

    std::size_t nInputGates = 0U;
    for (const auto& gate : *qcOriginal) {
        nInputGates++;
        mapGate(*gate);
        if (measureFrequency > 0 && nInputGates % measureFrequency == 0) {
            measureAndCorrect();
        }
    }

    // mapGate(...) can change 'isDecoded', therefore check it again
    if (!isDecoded) {
        measureAndCorrect();
        writeDecoding();
        isDecoded = true;
    }

    return qcMapped;
}
} // namespace ecc
