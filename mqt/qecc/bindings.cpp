/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md for more information.
 */

#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "DecodingSimulator.hpp"
#include "UFDecoder.hpp"
#include "UFHeuristic.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"

#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

std::vector<bool> sampleIidPauliErr(const std::size_t length, const double physicalErrRate) {
    return Utils::sampleErrorIidPauliNoise(length, physicalErrRate);
}

PYBIND11_MODULE(pyqecc, m) {
    m.doc() = "pybind11 for the MQT QECC quantum error-correcting codes tool";
    m.def("sample_iid_pauli_err", &sampleIidPauliErr, "Sample a iid pauli error represented as binary string");

    py::class_<Code>(m, "Code", "CSS code object")
            .def(py::init<>())
            .def(py::init<std::string, std::string>())
            .def(py::init<std::string>())
            .def(py::init<std::vector<std::vector<bool>>&>())
            .def(py::init<std::vector<std::vector<bool>>&, std::vector<std::vector<bool>>&>())
            .def("setHz", &Code::setHz)
            .def("getHz", &Code::getHzMat)
            .def("gethX", &Code::getHxMat)
            .def("getHx", &Code::setHx)
            .def_readwrite("n", &Code::n)
            .def_readwrite("k", &Code::k)
            .def_readwrite("d", &Code::d)
            .def("json", &Code::to_json)
            .def("is_x_stabilizer", &Code::isXStabilizer)
            .def("is_stabilizer", static_cast<bool (Code::*)(const std::vector<bool>&, const std::vector<bool>&) const>(&Code::isStabilizer))
            .def("is_stabilizer", static_cast<bool (Code::*)(const std::vector<bool>&) const>(&Code::isStabilizer))
            .def("get_syndrome", &Code::getSyndrome)
            .def("get_x_syndrome", &Code::getXSyndrome)
            .def("__repr__", &Code::toString);

    py::enum_<GrowthVariant>(m, "GrowthVariant")
            .value("ALL_COMPONENTS", GrowthVariant::AllComponents)
            .value("INVALID_COMPONENTS", GrowthVariant::InvalidComponents)
            .value("SINGLE_SMALLEST", GrowthVariant::SingleSmallest)
            .value("SINGLE_RANDOM", GrowthVariant::SingleRandom)
            .value("SINGLE_QUBIT_RANDOM", GrowthVariant::SingleQubitRandom)
            .export_values()
            .def(py::init([](const std::string& str) -> GrowthVariant { return growthVariantFromString(str); }));

    py::class_<DecodingResult>(m, "DecodingResult", "A decoding run result object")
            .def(py::init<>())
            .def_readwrite("decoding_time", &DecodingResult::decodingTime)
            .def_readwrite("estim_vec_idxs", &DecodingResult::estimNodeIdxVector)
            .def_readwrite("estimate", &DecodingResult::estimBoolVector)
            .def("json", &DecodingResult::to_json)
            .def("__repr__", &DecodingResult::toString);

    py::class_<Decoder>(m, "Decoder", "Decoder object")
            .def(py::init<>())
            .def_readwrite("result", &Decoder::result)
            .def_readwrite("growth", &Decoder::growth)
            .def("set_code", &Decoder::setCode)
            .def("set_growth", &Decoder::setGrowth)
            .def("decode", &Decoder::decode);

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

    py::class_<DecodingRunInformation>(m, "DecodingRunInformation")
            .def(py::init<>())
            .def(py::init<double, std::size_t, std::vector<bool>, std::vector<bool>, DecodingResultStatus, DecodingResult>())
            .def(py::init<double, std::size_t, std::vector<bool>, std::vector<bool>, DecodingResult>())
            .def_readwrite("physicalErrR", &DecodingRunInformation::physicalErrR)
            .def_readwrite("codeSize", &DecodingRunInformation::codeSize)
            .def_readwrite("error", &DecodingRunInformation::error)
            .def_readwrite("syndrome", &DecodingRunInformation::syndrome)
            .def_readwrite("status", &DecodingRunInformation::status)
            .def_readwrite("result", &DecodingRunInformation::result)
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

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
