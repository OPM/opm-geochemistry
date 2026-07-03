set(LIB_HEADERS GeoChemTestCase.hpp RegisterTests.hpp SetupAndRunSolver.hpp)
set(LIB_SOURCES RegisterTests.cpp SetupAndRunSolver.cpp)

add_library(RunGeoChem ${LIB_HEADERS} ${TEST_HEADERS} ${LIB_SOURCES})

target_link_libraries(RunGeoChem PRIVATE ${GEOCHEM_TEST_LIBRARY_TARGET})

target_include_directories(RunGeoChem
                           PRIVATE
                           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/opm/simulators/geochemistry/;${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                           )
