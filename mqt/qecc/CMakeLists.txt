# add pybind11 library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/pybind11" "extern/pybind11" EXCLUDE_FROM_ALL)

# add pybind11_json library
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/pybind11_json" "extern/pybind11_json" EXCLUDE_FROM_ALL)

pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib pybind11_json)

# LTO is disabled since a bug in manylinux2014 prevents the project from linking with it
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
    set_target_properties(${PROJECT_NAME}_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
endif ()
