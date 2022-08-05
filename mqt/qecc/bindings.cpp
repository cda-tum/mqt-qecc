/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md for more information.
 */

#include "Decoder.hpp"
#include "DecodingRunInformation.hpp"
#include "ImprovedUFD.hpp"
#include "OriginalUFD.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"
#include <pybind11/stl.h>

namespace py = pybind11;
namespace nl = nlohmann;
using namespace pybind11::literals;

std::vector<bool> sample_iid_pauli_err(const std::size_t length, const double physicalErrRate){
    return Utils::sampleErrorIidPauliNoise(length, physicalErrRate);
}

PYBIND11_MODULE(pyqecc, m) {
    m.doc() = "pybind11 for the MQT QECC quantum error-correcting codes tool";
    m.def("sample_iid_pauli_err", &sample_iid_pauli_err, "sample a iid pauli error represented as binary string");

    py::class_<Code>(m, "Code", "CSS code object")
            .def(py::init<>())
            .def(py::init<std::string>())
            .def(py::init<std::vector<std::vector<bool>>&>())
            .def("setHz", &Code::setHz)
            .def("getHz", &Code::getHzMat)
            .def_readwrite("N",&Code::N)
            .def_readwrite("K", &Code::K)
            .def_readwrite("D", &Code::D)
            .def("json", &Code::to_json)
            .def("get_syndrome", &Code::getSyndrome)
            .def("__repr__", &Code::toString);

    py::enum_<GrowthVariant>(m, "GrowthVariant")
            .value("ALL_COMPONENTS", GrowthVariant::ALL_COMPONENTS)
            .value("INVALID_COMPONENTS", GrowthVariant::INVALID_COMPONENTS)
            .value("SINGLE_SMALLEST", GrowthVariant::SINGLE_SMALLEST)
            .value("SINGLE_RANDOM", GrowthVariant::SINGLE_RANDOM)
            .value("SINGLE_QUBIT_RANDOM", GrowthVariant::SINGLE_QUBIT_RANDOM)
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

    py::class_<ImprovedUFD, Decoder>(m, "ImprovedUFD", "ImprovedUFD object")
            .def(py::init<>())
            .def_readwrite("result", &ImprovedUFD::result)
            .def_readwrite("growth", &ImprovedUFD::growth)
            .def("decode", &ImprovedUFD::decode);

    py::class_<OriginalUFD, Decoder>(m, "OriginalUFD", "OriginalUFD object")
            .def(py::init<>())
            .def_readwrite("result", &OriginalUFD::result)
            .def_readwrite("growth", &OriginalUFD::growth)
            .def("decode", &OriginalUFD::decode);

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
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}