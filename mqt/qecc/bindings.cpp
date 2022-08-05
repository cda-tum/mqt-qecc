/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md for more information.
 */

#include "Decoder.hpp"
#include "DecodingSimulator.hpp"
#include "ImprovedUFD.hpp"
#include "OriginalUFD.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"

namespace py = pybind11;
namespace nl = nlohmann;
using namespace pybind11::literals;

Code import_code_from_file(const std::string& filepath){
    return Code(filepath);
}

PYBIND11_MODULE(pyqecc, m) {
    m.doc() = "pybind11 for the MQT QECC quantum error-correcting codes tool";
    m.def("import_code_from_file", &import_code_from_file, "import a css code object from a file containing the parity-check matrix");

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
            .def("__repr__", &Code::toString);

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

}