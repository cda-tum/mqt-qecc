pybind11_add_module(
  pyqecc
  # Prefer thin LTO if available
  THIN_LTO
  # Optimize the bindings for size
  OPT_SIZE
  # Source code goes here
  bindings.cpp)
target_compile_definitions(pyqecc PRIVATE Z3_FOUND)
target_link_libraries(pyqecc PRIVATE MQT::QECC MQT::QECCFramework MQT::CorePython pybind11_json
                                     MQT::ProjectOptions MQT::ProjectWarnings)

# Install directive for scikit-build-core
install(TARGETS pyqecc LIBRARY DESTINATION . COMPONENT mqt-qecc_Python)
