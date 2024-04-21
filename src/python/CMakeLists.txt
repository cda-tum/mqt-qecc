target_compile_definitions(py${PROJECT_NAME} PRIVATE Z3_FOUND)

pybind11_add_module(
        pyqcec
        # Prefer thin LTO if available
        THIN_LTO
        # Optimize the bindings for size
        OPT_SIZE
        # Source code goes here
        bindings.cpp)
target_link_libraries(pyqcec PRIVATE MQT::QCEC MQT::CorePython pybind11_json
        MQT::ProjectOptions MQT::ProjectWarnings)

# Install directive for scikit-build-core
install(TARGETS pyqcec LIBRARY DESTINATION .)