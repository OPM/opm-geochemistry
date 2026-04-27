add_subdirectory(reference_results_eqsolver)
add_subdirectory(RunGeoChem)

set (CATCH2_SRC_DIR "${PROJECT_SOURCE_DIR}/tests/tests_using_catch2")
set (CATCH2_BIN_RUN_DIR "${CMAKE_BINARY_DIR}/bin/catch2_run_dir/")

set (EQSOLVER_REFERENCE_DIR "${CATCH2_SRC_DIR}/reference_results_eqsolver/refdata")

# Copy all reference simulation (archive) files to run directory
file(GLOB ALL_REFERENCE_CASES
     "${EQSOLVER_REFERENCE_DIR}/*.cereal"
)

file(MAKE_DIRECTORY "${CATCH2_BIN_RUN_DIR}")

foreach(CURRENT_FILE ${ALL_REFERENCE_CASES})
  get_filename_component(CURRENT_FILE_NAME ${CURRENT_FILE} NAME)
  configure_file("${CURRENT_FILE}" "${CATCH2_BIN_RUN_DIR}/${CURRENT_FILE_NAME}" COPYONLY)
endforeach()

set(CATCH2_SOURCE_FILES
    test_brent.cpp
    test_cubicrootsolver.cpp
    test_equilibrium_solver.cpp
    test_gauss_legendre.cpp
    test_geochem_phases.cpp
    test_input_reader.cpp
    test_input_json.cpp
    test_integration_for_diffuse_layer.cpp
    test_linear_algebra.cpp
    test_parse_chemistry.cpp
    test_string_functions.cpp
    test_thermo.cpp
    test_utility_functions.cpp
)

# Much fewer tests that are slower than the rest... only use when needed
set(CATCH2_SLOW_TESTS
    test_approval_transport.cpp
)

add_executable(catch2_tests catch2_main.cpp ${CATCH2_SOURCE_FILES} ${CATCH2_SLOW_TESTS})

if(Catch2_VERSION VERSION_GREATER_EQUAL 3.0)
  target_link_libraries(catch2_tests
   PRIVATE
     Catch2::Catch2WithMain
  )
else()
  target_link_libraries(catch2_tests
   PRIVATE
     Catch2::Catch2
  )
endif()

target_link_libraries(
    catch2_tests
    PRIVATE
    ApprovalTests::ApprovalTests
    ${GEOCHEM_TEST_LIBRARY_TARGET}
    EquilibriumTestCases
    RunGeoChem
)

set_target_properties(catch2_tests PROPERTIES CXX_STANDARD 17)

target_include_directories(catch2_tests
                           PRIVATE
                           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/opm/simulators/geochemistry/;${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                           )

add_test(NAME run_catch2_tests COMMAND catch2_tests WORKING_DIRECTORY ${CATCH2_BIN_RUN_DIR})

find_path(GEOCHEM_CATCH_MODULE_DIR Catch.cmake
          HINTS
            "${catch2_SOURCE_DIR}/contrib"
            "${Catch2_SOURCE_DIR}/contrib"
            "${catch2_SOURCE_DIR}/extras"
            "${Catch2_SOURCE_DIR}/extras")
if(GEOCHEM_CATCH_MODULE_DIR)
  list(APPEND CMAKE_MODULE_PATH ${GEOCHEM_CATCH_MODULE_DIR})
endif()
include(CTest)  # this adds BUILD_TESTING, which defaults to ON
include(Catch)
catch_discover_tests(catch2_tests WORKING_DIRECTORY ${CATCH2_BIN_RUN_DIR})
