add_subdirectory(EquilibriumTestCases)

if(GENERATE_NEW_REFERENCE_RESULTS)

    add_executable(GenRefResultsEQSOLVER generate_reference_results.cpp)

    target_link_libraries(
        GenRefResultsEQSOLVER
        PRIVATE
        ${GEOCHEM_TEST_LIBRARY_TARGET}
        EquilibriumTestCases
    )

    target_include_directories(GenRefResultsEQSOLVER
                               PRIVATE
                               "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>"
                               "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )

endif()
