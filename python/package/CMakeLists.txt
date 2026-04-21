# Create setup.py in build directory
configure_file(setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

add_custom_target(pygeochem-setuptools ALL
    DEPENDS project4py
    COMMAND ${CMAKE_COMMAND} -E rm -rf build  # remove build dir created by previous `python setup.py install` commands (see next) to ensure fresh rebuild since changed python files are not overwritten even with --force option
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pygeochem pygeochem
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:project4py> pygeochem
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${GEOCHEM_PYTHON_LIBRARY_TARGET}> pygeochem
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH="${PROJECT_PYTHONPATH}" ${CMAKE_COMMAND} -E true
    COMMAND ${PYTHON_EXECUTABLE} setup.py build --quiet --force
    COMMAND ${CMAKE_COMMAND} -E rm pygeochem/$<TARGET_FILE_NAME:project4py>
    COMMAND ${CMAKE_COMMAND} -E rm pygeochem/$<TARGET_FILE_NAME:${GEOCHEM_PYTHON_LIBRARY_TARGET}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# TODO: Not working it seems..
# Ensure that the path where the Python package is installed is non-empty
#if(NOT DEFINED PYTHON_INSTALL_PREFIX)
#    file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX} PYTHON_INSTALL_PREFIX)
#endif()

message("PYTHON_INSTALL_PREFIX ${PYTHON_INSTALL_PREFIX}")
message("{PROJECT_PYTHONPATH ${PROJECT_PYTHONPATH}")
message("PYTHON_SITE_PACKAGES_DIR ${PYTHON_SITE_PACKAGES_DIR}")

# Create an install target for the Python package
# NOTE: At least for now, remove the Python install prefix..
# (should work provided we have already activated a conda environment!)
install(CODE "execute_process(COMMAND pip install . --force-reinstall --no-compile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
# otherwise, installs into /usr/local by default...
#install(CODE "execute_process(COMMAND pip install . --no-deps --force-reinstall --no-compile --prefix=${PYTHON_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
