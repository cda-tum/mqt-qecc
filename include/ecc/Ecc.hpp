#pragma once

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ecc {
using Qubit      = qc::Qubit;
using QubitCount = std::size_t;

class Ecc {
public:
    enum class ID : std::uint8_t {
        Id,
        Q3Shor,
        Q9Shor,
        Q5Laflamme,
        Q7Steane,
        Q9Surface,
        Q18Surface
    };
    struct Info {
        ID          id;
        std::size_t nRedundantQubits; // usually number of physical qubits per
                                      // (encoded) logical qubit
        std::size_t nCorrectingBits;  // usually number of classical bits needed for
                                      // correcting one qubit
        std::string                                      name;
        std::vector<std::pair<std::size_t, std::string>> classicalRegisters;
    };

    Ecc(Info newEcc, std::shared_ptr<qc::QuantumComputation> qc,
        std::size_t newMeasureFrequency)
        : qcOriginal(std::move(qc)), measureFrequency(newMeasureFrequency),
          ecc(std::move(newEcc)) {
        qcMapped = std::make_shared<qc::QuantumComputation>();
    }
    virtual ~Ecc() = default;

    std::shared_ptr<qc::QuantumComputation> apply();

    [[nodiscard]] std::shared_ptr<qc::QuantumComputation>
    getOriginalCircuit() const {
        return qcOriginal;
    }

    [[nodiscard]] std::shared_ptr<qc::QuantumComputation>
    getMappedCircuit() const {
        return qcMapped;
    }

    virtual std::string getName() { return ecc.name; }

    virtual std::size_t getNOutputQubits(std::size_t nInputQubits) {
        return nInputQubits * ecc.nRedundantQubits + ecc.nCorrectingBits;
    }

    [[nodiscard]] std::vector<Qubit> getDataQubits() const {
        auto               numberOfDataQubits = qcOriginal->getNqubits() * ecc.nRedundantQubits;
        std::vector<Qubit> dataQubits(numberOfDataQubits);
        std::iota(std::begin(dataQubits), std::end(dataQubits), 0);
        return dataQubits;
    }

    std::shared_ptr<qc::QuantumComputation> qcOriginal;
    std::shared_ptr<qc::QuantumComputation> qcMapped;
    std::size_t                             measureFrequency;
    bool                                    isDecoded    = true;
    bool                                    gatesWritten = false;
    Info                                    ecc;

protected:
    void initMappedCircuit();

    /**
     * prepares an encoded logical |0> state in the qcMapped circuit.
     * May, but does not have to be overridden by subclasses.
     * */
    virtual void writeEncoding() {
        if (!isDecoded) {
            return;
        }
        isDecoded = false;
        measureAndCorrect();
    }

    /**
     * in case of an error, calling this function creates a 'clean' state again.
     * Usual structure:
     *
     * for each logical qubit i:
     * -- reset ancilla qubits
     * -- measure physical data qubits of logical qubit[i] onto ancilla qubits
     * -- correct data qubits based on measurement results of ancilla qubits
     * */
    virtual void measureAndCorrect() = 0;

    /**
     * moves encoded state information back to original qubits (i.e. 1 qubit per
     * logical qubit)
     * */
    virtual void writeDecoding() = 0;

    virtual void mapGate(const qc::Operation& gate) = 0;

    void gateNotAvailableError(const qc::Operation& gate) const {
        std::stringstream stream;
        stream << "Gate " << toString(gate.getType())
               << " not supported to encode in error code " << ecc.name << "!";
        throw qc::QFRException(stream.str());
    }
};
} // namespace ecc
