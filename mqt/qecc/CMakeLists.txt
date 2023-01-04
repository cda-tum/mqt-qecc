# add MQT::QFR library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr" EXCLUDE_FROM_ALL)

pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib pybind11_json MQT::qfr_ecc)
target_compile_definitions(py${PROJECT_NAME} PRIVATE VERSION_INFO=${QECC_VERSION_INFO})

target_compile_definitions(py${PROJECT_NAME} PUBLIC VERSION_INFO=${QECC_VERSION_INFO})
