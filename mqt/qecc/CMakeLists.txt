pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_qecc_lib MQT::qfr_python pybind11_json)

# LTO is disabled since a bug in manylinux2014 prevents the project from linking with it
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
    set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
    set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
    set_target_properties(qfr PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
    set_target_properties(DDPackage PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
endif ()
