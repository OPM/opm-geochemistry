file(GLOB_RECURSE HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
file(GLOB_RECURSE CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_library(EquilibriumTestCases ${HPP_FILES} ${CPP_FILES})

target_link_libraries(EquilibriumTestCases PRIVATE ${GEOCHEM_TEST_LIBRARY_TARGET})

target_include_directories(EquilibriumTestCases
                           PRIVATE
                           "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/opm/simulators/geochemistry/;${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                           )
