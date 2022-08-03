/*
 * This file is part of MQT QECC library which is released under the MIT license.
 * See file README.md or go to todo for more information.
 */

#include "Decoder.hpp"
#include "nlohmann/json.hpp"
#include "pybind11/pybind11.h"
#include "pybind11_json/pybind11_json.hpp"

namespace py = pybind11;
namespace nl = nlohmann;
using namespace pybind11::literals;
//todo

Code import_code_from_file(const std::string& filepath){
    return Code(filepath);
}

PYBIND11_MODULE(pyqecc, m) {
    m.doc() = "pybind11 for the MQT QECC quantum error-correcting codes tool";
    m.def("import_code_from_file", &import_code_from_file, "import a css code object from a file containing the parity-check matrix");

    py::class_<Code>(m, "Code", "CSS code object")
            .def(py::init<>())
            .def_readwrite("Hz", &Code::Hz)
            .def_readwrite("N",&Code::N)
            .def_readwrite("K", &Code::K)
            .def_readwrite("D", &Code::D)
            .def("json", &Code::to_json)
            .def("__repr__", &Code::toString)
}