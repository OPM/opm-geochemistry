if(NOT GEOCHEM_PYTHON_LIBRARY_TARGET)
    message(FATAL_ERROR "GEOCHEM_PYTHON_LIBRARY_TARGET must be set before adding the python subdirectory.")
endif()

add_subdirectory(package)

# The Python extension links this library into a shared module.
set_target_properties(${GEOCHEM_PYTHON_LIBRARY_TARGET}
    PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Collect all .py.cxx files, and add them to Python module
# file(GLOB_RECURSE CXX_FILES_PY RELATIVE ${PROJECT_SOURCE_DIR}/opm/simulators/ *.py.cxx)
# pybind11_add_module(project4py ${CXX_FILES_PY})
pybind11_add_module(project4py ${PROJECT_SOURCE_DIR}/opm/simulators/GeochemPythonInterface.py.cxx)

# Add the root directory of the project to the include list
target_include_directories(project4py PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/opm/simulators/geochemistry/)

# Link libraries
target_link_libraries(project4py PRIVATE ${GEOCHEM_PYTHON_LIBRARY_TARGET})

# This is needed to avoid compilation issues when using clang-9 and c++17
if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0.0)
    target_compile_options(project4py PRIVATE -fsized-deallocation)
endif()
