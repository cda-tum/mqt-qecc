pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib pybind11_json MQT::qfr_ecc
                                               MQT::qfr_python)
target_compile_definitions(py${PROJECT_NAME} PRIVATE VERSION_INFO=${QECC_VERSION_INFO})
