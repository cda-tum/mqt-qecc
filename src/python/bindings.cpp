/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md for more information.
 */

#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "UFDecoder.hpp"
#include "UFHeuristic.hpp"
#include "ecc/Ecc.hpp"
#include "ecc/Id.hpp"
#include "ecc/Q18Surface.hpp"
#include "ecc/Q3Shor.hpp"
#include "ecc/Q5Laflamme.hpp"
#include "ecc/Q7Steane.hpp"
#include "ecc/Q9Shor.hpp"
#include "ecc/Q9Surface.hpp"
#include "ir/QuantumComputation.hpp"
#include "python/qiskit/QuantumCircuit.hpp"

#include <pybind11/pybind11.h>
#include <pybind11_json/pybind11_json.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

std::vector<bool> sampleIidPauliErr(const std::size_t length, const double physicalErrRate) {
    return Utils::sampleErrorIidPauliNoise(length, physicalErrRate);
}

py::dict applyEcc(const py::object& circ, const std::string& eccName, const size_t eccFrequency) {
    auto qc = std::make_shared<qc::QuantumComputation>();

    if (py::isinstance<py::str>(circ)) {
        auto&& file = circ.cast<std::string>();
        qc->import(file);
    } else {
        py::object const quantumCircuit = py::module::import("qiskit").attr("QuantumCircuit");
        if (py::isinstance(circ, quantumCircuit)) {
            qc::qiskit::QuantumCircuit::import(*qc, circ);
        }
    }

    std::unique_ptr<ecc::Ecc> mapper{};

    if (eccName == "Id") {
        mapper = std::make_unique<ecc::Id>(qc, eccFrequency);
    } else if (eccName == "Q3Shor") {
        mapper = std::make_unique<ecc::Q3Shor>(qc, eccFrequency);
    } else if ((eccName == "Q5Laflamme")) {
        mapper = std::make_unique<ecc::Q5Laflamme>(qc, eccFrequency);
    } else if ((eccName == "Q7Steane")) {
        mapper = std::make_unique<ecc::Q7Steane>(qc, eccFrequency);
    } else if ((eccName == "Q9Shor")) {
        mapper = std::make_unique<ecc::Q9Shor>(qc, eccFrequency);
    } else if ((eccName == "Q9Surface")) {
        mapper = std::make_unique<ecc::Q9Surface>(qc, eccFrequency);
    } else if ((eccName == "Q18Surface")) {
        mapper = std::make_unique<ecc::Q18Surface>(qc, eccFrequency);
    } else {
        std::stringstream ss{};
        ss << "No ECC found for " << eccName << " ";
        ss << "Available ECCs: ";
        ss << "Id"
           << ", ";
        ss << "Q3Shor"
           << ", ";
        ss << "Q5Laflamme"
           << ", ";
        ss << "Q7Steane"
           << ", ";
        ss << "Q9Shor"
           << ", ";
        ss << "Q9Surface"
           << ", ";
        ss << "Q18Surface";
        throw std::invalid_argument(ss.str());
    }

    std::shared_ptr<qc::QuantumComputation> const qcECC = mapper->apply();

    std::ostringstream oss{};
    qcECC->dump(oss, qc::Format::OpenQASM2);

    return py::dict("circ"_a = oss.str());
}

PYBIND11_MODULE(pyqecc, m, py::mod_gil_not_used()) {
    // In this function the exposed python functions are defined.
    // Additionally, the required parameters are described.
    m.doc() = "pybind11 for the MQT QECC quantum error-correcting codes tool";
    m.def("sample_iid_pauli_err", &sampleIidPauliErr, "Sample a iid pauli error represented as binary string");

    py::class_<Code>(m, "Code", "CSS code object")
            .def(py::init<>(), "Constructs Code object without pcms set")
            .def(py::init<std::string, std::string>(), "Constructs Code object with pcms read from file paths")
            .def(py::init<std::string>(), "Constructs single-sided CSS Code object with pcm read from file path")
            .def(py::init<std::vector<std::vector<bool>>&>(), "Constructs single-sided CSS Code object with pcm given as boolean matrix")
            .def(py::init<std::vector<std::vector<bool>>&, std::vector<std::vector<bool>>&>(), "Constructs single-sided CSS Code object with pcms given as boolean matrices")
            .def("setHz", &Code::setHz)
            .def("getHz", &Code::getHzMat)
            .def("setHx", &Code::setHx)
            .def("getHx", &Code::getHxMat)
            .def_readwrite("n", &Code::n)
            .def_readwrite("k", &Code::k)
            .def_readwrite("d", &Code::d)
            .def("json", &Code::to_json)
            .def("is_x_stabilizer", &Code::isXStabilizer, "Checks if vector is a x-stabilizer")
            .def("is_stabilizer", static_cast<bool (Code::*)(const std::vector<bool>&, const std::vector<bool>&) const>(&Code::isStabilizer), "Checks if vector is a stabilizer")
            .def("is_stabilizer", static_cast<bool (Code::*)(const std::vector<bool>&) const>(&Code::isStabilizer))
            .def("get_syndrome", &Code::getSyndrome, "Computes syndrome vector")
            .def("get_x_syndrome", &Code::getXSyndrome, "Computes single sided syndrome vector (length n)")
            .def("__repr__", &Code::toString);

    py::enum_<GrowthVariant>(m, "GrowthVariant")
            .value("ALL_COMPONENTS", GrowthVariant::AllComponents, "Grows all components in each iteration")
            .value("INVALID_COMPONENTS", GrowthVariant::InvalidComponents, "Grows invalid components only in each iteration")
            .value("SINGLE_SMALLEST", GrowthVariant::SingleSmallest, "Grows only smallest component in each iteration")
            .value("SINGLE_RANDOM", GrowthVariant::SingleRandom, "Grows a single uniformly random component in each iteration")
            .value("SINGLE_QUBIT_RANDOM", GrowthVariant::SingleQubitRandom, "Grows component around a single qubit in each iteration")
            .export_values()
            .def(py::init([](const std::string& str) -> GrowthVariant { return growthVariantFromString(str); }));

    py::class_<DecodingResult>(m, "DecodingResult", "Holds information about a single decoding step")
            .def(py::init<>())
            .def_readwrite("decoding_time", &DecodingResult::decodingTime, "Time used for a single decoding step")
            .def_readwrite("estim_vec_idxs", &DecodingResult::estimNodeIdxVector, "Computed estimates given as indices (over qubits)")
            .def_readwrite("estimate", &DecodingResult::estimBoolVector, "Computed estimate as boolean vector")
            .def("json", &DecodingResult::to_json)
            .def("__repr__", &DecodingResult::toString);

    py::class_<Decoder>(m, "Decoder", "Decoder object")
            .def(py::init<>())
            .def_readwrite("result", &Decoder::result, "Decoding result object")
            .def_readwrite("growth", &Decoder::growth, "The growth variant currently set")
            .def("set_code", &Decoder::setCode)
            .def("set_growth", &Decoder::setGrowth)
            .def("decode", &Decoder::decode, "Decode a syndrome vector. After completion the result field is not null");

    py::class_<UFHeuristic, Decoder>(m, "UFHeuristic", "UFHeuristic object")
            .def(py::init<>())
            .def_readwrite("result", &UFHeuristic::result)
            .def_readwrite("growth", &UFHeuristic::growth)
            .def("reset", &UFHeuristic::reset)
            .def("decode", &UFHeuristic::decode);

    py::class_<UFDecoder, Decoder>(m, "UFDecoder", "UFDecoder object")
            .def(py::init<>())
            .def_readwrite("result", &UFDecoder::result)
            .def_readwrite("growth", &UFDecoder::growth)
            .def("decode", &UFDecoder::decode);

    py::enum_<DecodingResultStatus>(m, "DecodingResultStatus")
            .value("ALL_COMPONENTS", DecodingResultStatus::SUCCESS)
            .value("INVALID_COMPONENTS", DecodingResultStatus::FAILURE)
            .export_values()
            .def(py::init([](const std::string& str) -> DecodingResultStatus { return decodingResultStatusFromString(str); }));

    py::class_<DecodingRunInformation>(m, "DecodingRunInformation", "Holds additional information about a decoding run")
            .def(py::init<>())
            .def(py::init<double, std::size_t, std::vector<bool>, std::vector<bool>, DecodingResultStatus, DecodingResult>())
            .def(py::init<double, std::size_t, std::vector<bool>, std::vector<bool>, DecodingResult>())
            .def_readwrite("physicalErrR", &DecodingRunInformation::physicalErrR, "Physical error rate")
            .def_readwrite("codeSize", &DecodingRunInformation::codeSize, "Code size")
            .def_readwrite("error", &DecodingRunInformation::error, "Error vector")
            .def_readwrite("syndrome", &DecodingRunInformation::syndrome, "Syndrome vector")
            .def_readwrite("status", &DecodingRunInformation::status, "Result status (success or failure)")
            .def_readwrite("result", &DecodingRunInformation::result, "Result object")
            .def("json", &DecodingRunInformation::to_json)
            .def("__repr__", &DecodingRunInformation::toString);

    py::class_<DecodingSimulator>(m, "DecodingSimulator")
            .def(py::init<>())
            .def("simulate_wer", &DecodingSimulator::simulateWER)
            .def("simulate_avg_runtime", &DecodingSimulator::simulateAverageRuntime);

    py::enum_<DecoderType>(m, "DecoderType")
            .value("UF_HEURISTIC", DecoderType::UfHeuristic)
            .value("ORIGINAL_UF", DecoderType::UfDecoder)
            .export_values()
            .def(py::init([](const std::string& str) -> DecoderType { return decoderTypeFromString(str); }));

    // Function signature for applying an error correction scheme
    m.def("apply_ecc", &applyEcc, "apply an ECC to a circuit and return the resulting circuit as an OpenQASM string",
          "circuit_name"_a,
          "ecc_name"_a,
          "ecc_frequency"_a = 100);
}
