# add pybind11 library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/pybind11" "extern/pybind11" EXCLUDE_FROM_ALL)

# add pybind11_json library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/pybind11_json" "extern/pybind11_json" EXCLUDE_FROM_ALL)

pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib pybind11_json)
target_compile_definitions(py${PROJECT_NAME} PRIVATE VERSION_INFO=${QECC_VERSION_INFO})
