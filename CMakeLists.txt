###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 3.23)
project(opm-geochemistry C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

include(GNUInstallDirs)
include(CMakeLists_files.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Other options
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(BUILD_PYTHON "Attempt to build Python package?" OFF)
if(WIN32)
  set(_BUILD_STANDALONE_DEFAULT ON)
  set(_ENABLE_OPM_INTEGRATION_DEFAULT OFF)
else()
  set(_BUILD_STANDALONE_DEFAULT OFF)
  set(_ENABLE_OPM_INTEGRATION_DEFAULT ON)
endif()
option(BUILD_STANDALONE "Build the standalone geochemistry library and tools?" ${_BUILD_STANDALONE_DEFAULT})
option(ENABLE_OPM_INTEGRATION "Build the OPM integration layer?" ${_ENABLE_OPM_INTEGRATION_DEFAULT})

if(WIN32 AND ENABLE_OPM_INTEGRATION)
  message(WARNING "OPM integration is usually unavailable on Windows. Configure with -DENABLE_OPM_INTEGRATION=OFF for standalone builds.")
endif()

if(NOT BUILD_STANDALONE AND NOT ENABLE_OPM_INTEGRATION)
  message(FATAL_ERROR "Nothing to build. Enable BUILD_STANDALONE and/or ENABLE_OPM_INTEGRATION.")
endif()

function(geochem_setup_python_dependencies)
  if(PYTHON_EXECUTABLE AND NOT Python3_EXECUTABLE)
    set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
  endif()

  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  if(Python3_VERSION_MINOR LESS 3)
    message(SEND_ERROR "OPM requires python >= 3.3 but only version ${Python3_VERSION} was found")
  endif()

  if(NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH "Python interpreter" FORCE)
  endif()

  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    include(DownloadPyBind11)
  endif()
endfunction()

function(geochem_setup_python_packaging)
  execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
          from distutils import sysconfig as sc
          print(sc.get_python_lib(prefix='', plat_specific=True))"
      OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)

  file(TO_NATIVE_PATH ${CMAKE_BINARY_DIR}/python/package/build/lib PYTHON_PACKAGE_DIR_1)
  file(TO_NATIVE_PATH ${CMAKE_BINARY_DIR}/python/package PYTHON_PACKAGE_DIR_2)

  if(WIN32)
    set(PROJECT_PYTHONPATH "${PYTHON_PACKAGE_DIR_1}\;${PYTHON_PACKAGE_DIR_2}\;$ENV{PYTHONPATH}" PARENT_SCOPE)
  else()
    set(PROJECT_PYTHONPATH '${PYTHON_PACKAGE_DIR_1}:${PYTHON_PACKAGE_DIR_2}:$ENV{PYTHONPATH}' PARENT_SCOPE)
  endif()

  set(PYTHON_SITE_PACKAGES_DIR "${PYTHON_SITE_PACKAGES_DIR}" PARENT_SCOPE)
  set(PYTHON_PACKAGE_DIR_1 "${PYTHON_PACKAGE_DIR_1}" PARENT_SCOPE)
  set(PYTHON_PACKAGE_DIR_2 "${PYTHON_PACKAGE_DIR_2}" PARENT_SCOPE)
endfunction()

function(geochem_setup_test_dependencies)
  find_package(ApprovalTests CONFIG QUIET)
  if(NOT ApprovalTests_FOUND)
    if(EXISTS "${PROJECT_SOURCE_DIR}/build/_deps/approvaltests-src/CMakeLists.txt")
      add_subdirectory(
        "${PROJECT_SOURCE_DIR}/build/_deps/approvaltests-src"
        "${CMAKE_BINARY_DIR}/_deps/approvaltests-build"
        EXCLUDE_FROM_ALL
      )
    else()
      include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindApprovalTests.cmake)
    endif()
  endif()

  find_package(Catch2 CONFIG QUIET)
  if(NOT Catch2_FOUND)
    if(EXISTS "${PROJECT_SOURCE_DIR}/build/_deps/catch2-src/CMakeLists.txt")
      add_subdirectory(
        "${PROJECT_SOURCE_DIR}/build/_deps/catch2-src"
        "${CMAKE_BINARY_DIR}/_deps/catch2-build"
        EXCLUDE_FROM_ALL
      )
      set(catch2_SOURCE_DIR "${PROJECT_SOURCE_DIR}/build/_deps/catch2-src" PARENT_SCOPE)
    else()
      include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindCatch2.cmake)
    endif()
  endif()
endfunction()

function(geochem_ensure_standalone_dependencies)
  find_package(fmt CONFIG QUIET)
  if(TARGET fmt::fmt-header-only)
    set(GEOCHEM_FMT_TARGET fmt::fmt-header-only PARENT_SCOPE)
  elseif(TARGET fmt::fmt)
    set(GEOCHEM_FMT_TARGET fmt::fmt PARENT_SCOPE)
  else()
    find_path(GEOCHEM_FMT_INCLUDE_DIR fmt/core.h)
    if(GEOCHEM_FMT_INCLUDE_DIR)
      if(NOT TARGET geochem_fmt_headers)
        add_library(geochem_fmt_headers INTERFACE)
        target_include_directories(geochem_fmt_headers INTERFACE "${GEOCHEM_FMT_INCLUDE_DIR}")
        target_compile_definitions(geochem_fmt_headers INTERFACE FMT_HEADER_ONLY)
      endif()
      if(NOT TARGET fmt::fmt-header-only)
        add_library(fmt::fmt-header-only ALIAS geochem_fmt_headers)
      endif()
      set(GEOCHEM_FMT_TARGET fmt::fmt-header-only PARENT_SCOPE)
    else()
      include(Fetchfmt)
      if(TARGET fmt::fmt-header-only)
        set(GEOCHEM_FMT_TARGET fmt::fmt-header-only PARENT_SCOPE)
      elseif(TARGET fmt::fmt)
        set(GEOCHEM_FMT_TARGET fmt::fmt PARENT_SCOPE)
      else()
        message(FATAL_ERROR "Could not locate fmt.")
      endif()
    endif()
  endif()

  find_package(cereal CONFIG QUIET)
  if(TARGET cereal::cereal)
    set(GEOCHEM_CEREAL_TARGET cereal::cereal PARENT_SCOPE)
  else()
    find_path(GEOCHEM_CEREAL_INCLUDE_DIR cereal/cereal.hpp)
    if(NOT GEOCHEM_CEREAL_INCLUDE_DIR)
      include(Findcereal)
      set(GEOCHEM_CEREAL_INCLUDE_DIR "${cereal_SOURCE_DIR}/include")
    endif()
    if(NOT GEOCHEM_CEREAL_INCLUDE_DIR)
      message(FATAL_ERROR "Could not locate cereal.")
    endif()
    if(NOT TARGET geochem_cereal_headers)
      add_library(geochem_cereal_headers INTERFACE)
      target_include_directories(geochem_cereal_headers INTERFACE "${GEOCHEM_CEREAL_INCLUDE_DIR}")
    endif()
    if(NOT TARGET cereal::cereal)
      add_library(cereal::cereal ALIAS geochem_cereal_headers)
    endif()
    set(GEOCHEM_CEREAL_TARGET cereal::cereal PARENT_SCOPE)
  endif()

  find_package(nlohmann_json CONFIG QUIET)
  if(TARGET nlohmann_json::nlohmann_json)
    set(GEOCHEM_JSON_TARGET nlohmann_json::nlohmann_json PARENT_SCOPE)
  elseif(TARGET nlohmann_json)
    set(GEOCHEM_JSON_TARGET nlohmann_json PARENT_SCOPE)
  else()
    find_path(GEOCHEM_JSON_INCLUDE_DIR nlohmann/json.hpp
      HINTS
        "${PROJECT_SOURCE_DIR}/build/_deps/json-src/include"
        "${PROJECT_BINARY_DIR}/_deps/json-src/include"
    )
    if(GEOCHEM_JSON_INCLUDE_DIR)
      if(NOT TARGET geochem_nlohmann_json_headers)
        add_library(geochem_nlohmann_json_headers INTERFACE)
        target_include_directories(geochem_nlohmann_json_headers INTERFACE "${GEOCHEM_JSON_INCLUDE_DIR}")
      endif()
      if(NOT TARGET nlohmann_json::nlohmann_json)
        add_library(nlohmann_json::nlohmann_json ALIAS geochem_nlohmann_json_headers)
      endif()
      set(GEOCHEM_JSON_TARGET nlohmann_json::nlohmann_json PARENT_SCOPE)
    else()
      include(Findnlohmann_json)
      if(TARGET nlohmann_json::nlohmann_json)
        set(GEOCHEM_JSON_TARGET nlohmann_json::nlohmann_json PARENT_SCOPE)
      elseif(TARGET nlohmann_json)
        set(GEOCHEM_JSON_TARGET nlohmann_json PARENT_SCOPE)
      else()
        message(FATAL_ERROR "Could not locate nlohmann_json.")
      endif()
    endif()
  endif()
endfunction()

function(geochem_add_standalone_targets)
  geochem_ensure_standalone_dependencies()

  add_library(GeoChem
    ${GEOCHEM_CORE_PUBLIC_HEADER_FILES}
    ${GEOCHEM_CORE_SOURCE_FILES}
  )
  add_library(GeoChem::GeoChem ALIAS GeoChem)

  target_compile_features(GeoChem PUBLIC cxx_std_17)
  target_include_directories(GeoChem
    PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_link_libraries(GeoChem
    PUBLIC
      ${GEOCHEM_FMT_TARGET}
      ${GEOCHEM_CEREAL_TARGET}
      ${GEOCHEM_JSON_TARGET}
  )

  find_package(netCDFCxx CONFIG QUIET)
  if(netCDFCxx_FOUND)
    message(STATUS "Found NetCDF for standalone GeoChem.")
    target_link_libraries(GeoChem PUBLIC netCDF::netcdf netCDF::netcdf-cxx4)
    target_compile_definitions(GeoChem PUBLIC NETCDF=1)
  else()
    message(STATUS "Did not find NetCDF for standalone GeoChem, not using it...")
  endif()

  add_executable(GeoChemX examples/GeoChemX.cpp)
  target_link_libraries(GeoChemX PRIVATE GeoChem)

  add_executable(ThermoX examples/ThermoX.cpp)
  target_link_libraries(ThermoX PRIVATE GeoChem)
endfunction()

macro(opm-geochemistry_prereqs_hook)
  target_link_libraries(opmgeochemistry
    PUBLIC
      opmsimulators
      opmcommon
      opmgrid
  )

  find_package(ApprovalTests CONFIG)
  if (NOT ApprovalTests_FOUND)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindApprovalTests.cmake)
  endif()
  find_package(Catch2 CONFIG)
  if (NOT Catch2_FOUND)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindCatch2.cmake)
  endif()
  find_package(cereal CONFIG QUIET)
  if (NOT cereal_FOUND)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/Findcereal.cmake)
  endif()
  find_package(nlohmann_json CONFIG QUIET)
  if (NOT nlohmann_json_FOUND)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/Findnlohmann_json.cmake)
  endif()

  if (BUILD_PYTHON AND NOT BUILD_STANDALONE)
    geochem_setup_python_dependencies()
  endif()
endmacro()

macro(opm-geochemistry_tests_hook)
  if(BUILD_TESTING)
    set(GEOCHEM_TEST_LIBRARY_TARGET opmgeochemistry)
    add_subdirectory(tests)
  endif()
endmacro()

macro(opm-geochemistry_targets_hook)
  opm_add_library(
    TARGET
      moduleVersionGeochem
    TYPE
      OBJECT
    SOURCES
      opm/simulators/utils/moduleVersion.cpp
  )
  add_dependencies(moduleVersionGeochem opmgeochemistry)

  opm_add_test(flow_blackoil_geochemistry
               ONLY_COMPILE
               ALWAYS_ENABLE
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
               SOURCES
                 flow/flow_blackoil_geochemistry.cpp
                 $<TARGET_OBJECTS:moduleVersionGeochem>
  )

  opm_add_test(flow_onephase_geochemistry
               ONLY_COMPILE
               ALWAYS_ENABLE
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
               SOURCES
                 flow/flow_onephase_geochemistry.cpp
                 $<TARGET_OBJECTS:moduleVersionGeochem>
  )

  opm_add_test(test_opmgeocheminterface
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
                 Boost::unit_test_framework
               SOURCES
                 tests/test_opmgeocheminterface.cpp
                 $<TARGET_OBJECTS:moduleVersionGeochem>
  )

  if(NOT BUILD_STANDALONE)
    add_executable(GeoChemX
      examples/GeoChemX.cpp
    )
    target_link_libraries(GeoChemX
      PRIVATE
        opmgeochemistry
    )

    add_executable(ThermoX
      examples/ThermoX.cpp
    )
    target_link_libraries(ThermoX
      PRIVATE
        opmgeochemistry
    )
  endif()

  if(BUILD_PYTHON AND NOT BUILD_STANDALONE)
    set(GEOCHEM_PYTHON_LIBRARY_TARGET opmgeochemistry)
    geochem_setup_python_packaging()
    add_subdirectory(python)
  endif()
endmacro()

if(BUILD_STANDALONE)
  geochem_add_standalone_targets()
  set(GEOCHEM_PYTHON_LIBRARY_TARGET GeoChem)

  include(CTest)
  if(BUILD_TESTING)
    set(GEOCHEM_TEST_LIBRARY_TARGET GeoChem)
    geochem_setup_test_dependencies()
    add_subdirectory(tests)
  endif()

  if(BUILD_PYTHON)
    geochem_setup_python_dependencies()
    geochem_setup_python_packaging()
    add_subdirectory(python)
  endif()
endif()

if(ENABLE_OPM_INTEGRATION)
  if(SIBLING_SEARCH AND NOT opm-common_DIR)
    get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
    get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
    get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
    get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
    if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
      set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
      set(opm-grid_DIR ${_modules_dir}/opm-grid/${_leaf_dir_name})
      set(opm-simulators_DIR ${_modules_dir}/opm-simulators/${_leaf_dir_name})
    else()
      string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
      if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
          AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
        set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
        string(REPLACE ${PROJECT_NAME} opm-grid _opm_grid_leaf ${_leaf_dir_name})
        set(opm-grid_DIR ${_parent_full_dir}/${_opm_grid_leaf})
        string(REPLACE ${PROJECT_NAME} opm-simulators _opm_simulators_leaf ${_leaf_dir_name})
        set(opm-simulators_DIR ${_parent_full_dir}/${_opm_grid_leaf})
      elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
        set(opm-common_DIR "${_parent_full_dir}/opm-common")
        set(opm-grid_DIR "${_parent_full_dir}/opm-grid")
        set(opm-simulators_DIR "${_parent_full_dir}/opm-grid")
      endif()
    endif()
  endif()
  if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
    message(WARNING "Value ${opm-common_DIR} passed to variable"
      " opm-common_DIR is not a directory")
  endif()

  find_package(opm-common REQUIRED)

  include (OpmInit)
  include (OpmLibMain)
endif()
