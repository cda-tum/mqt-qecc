/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md for more information.
 */

#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "UFDecoder.hpp"
#include "UFHeuristic.hpp"
#include "eccs/Ecc.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"

#include <../extern/qfr/include/QuantumComputation.hpp>
#include <../extern/qfr/mqt/qfr/qiskit/QasmQobjExperiment.hpp>
#include <../extern/qfr/mqt/qfr/qiskit/QuantumCircuit.hpp>
#include <eccs/IdEcc.hpp>
#include <eccs/Q18SurfaceEcc.hpp>
#include <eccs/Q3ShorEcc.hpp>
#include <eccs/Q5LaflammeEcc.hpp>
#include <eccs/Q7SteaneEcc.hpp>
#include <eccs/Q9ShorEcc.hpp>
#include <eccs/Q9SurfaceEcc.hpp>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

std::vector<bool> sampleIidPauliErr(const std::size_t length, const double physicalErrRate) {
    return Utils::sampleErrorIidPauliNoise(length, physicalErrRate);
}

py::dict apply_ecc(const py::object& circ, const std::string& eccString, const int ecc_frequency) {
    qc::QuantumComputation qc{};
    std::string            eccName{eccString};

    try {
        if (py::isinstance<py::str>(circ)) {
            auto&& file = circ.cast<std::string>();
            qc.import(file);
        } else {
            py::object QuantumCircuit       = py::module::import("qiskit").attr("QuantumCircuit");
            py::object pyQasmQobjExperiment = py::module::import("qiskit.qobj").attr("QasmQobjExperiment");
            if (py::isinstance(circ, QuantumCircuit)) {
                qc::qiskit::QuantumCircuit::import(qc, circ);
            } else if (py::isinstance(circ, pyQasmQobjExperiment)) {
                qc::qiskit::QasmQobjExperiment::import(qc, circ);
            }
        }
    } catch (std::exception const& e) {
        std::stringstream ss{};
        ss << "Could not import circuit: " << e.what();
        return py::dict("error"_a = ss.str());
    }

    Ecc* mapper           = nullptr;
    int  measureFrequency = ecc_frequency;

    if (eccName.compare(IdEcc::getName()) == 0) {
        mapper = new IdEcc(qc, measureFrequency);
    } else if (eccName.compare(Q3ShorEcc::getName()) == 0) {
        mapper = new Q3ShorEcc(qc, measureFrequency);
    } else if (eccName.compare(Q5LaflammeEcc::getName()) == 0) {
        mapper = new Q5LaflammeEcc(qc, measureFrequency);
    } else if (eccName.compare(Q7SteaneEcc::getName()) == 0) {
        mapper = new Q7SteaneEcc(qc, measureFrequency);
    } else if (eccName.compare(Q9ShorEcc::getName()) == 0) {
        mapper = new Q9ShorEcc(qc, measureFrequency);
    } else if (eccName.compare(Q9SurfaceEcc::getName()) == 0) {
        mapper = new Q9SurfaceEcc(qc, measureFrequency);
    } else if (eccName.compare(Q18SurfaceEcc::getName()) == 0) {
        mapper = new Q18SurfaceEcc(qc, measureFrequency);
    } else {
        std::stringstream ss{};
        ss << "No ECC found for " << eccName << " ";
        ss << "Available ECCs: ";
        ss << IdEcc::getName() << ", ";
        ss << Q3ShorEcc::getName() << ", ";
        ss << Q5LaflammeEcc::getName() << ", ";
        ss << Q7SteaneEcc::getName() << ", ";
        ss << Q9ShorEcc::getName() << ", ";
        ss << Q9SurfaceEcc::getName() << ", ";
        ss << Q18SurfaceEcc::getName();
        return py::dict("error"_a = ss.str());
    }
    qc::QuantumComputation& qcECC = mapper->apply();

    std::ostringstream oss{};
    qcECC.dump(oss, qc::OpenQASM);

    return py::dict("circ"_a = oss.str());
}

PYBIND11_MODULE(pyqecc, m) {
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

    m.def("apply_ecc", &apply_ecc, "applying an ecc to a circuit an returning a openQasm dump",
          "circ"_a,
          "eccString"_a,
          "ecc_frequency"_a = 100);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
