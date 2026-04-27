# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file separates the standalone geochemistry solver from the optional
# OPM integration layer so both build paths can reuse the same source lists.

set(GEOCHEM_CORE_SOURCE_FILES
    opm/simulators/geochemistry/GeoChemIF.cpp
    opm/simulators/geochemistry/StandaloneSolvers.cpp
    opm/simulators/geochemistry/Common/ChemGlobal.cpp
    opm/simulators/geochemistry/Core/ChemBasVec.cpp
    opm/simulators/geochemistry/Core/ChemGCSolver.cpp
    opm/simulators/geochemistry/Core/ChemInitChem.cpp
    opm/simulators/geochemistry/Core/ChemTable.cpp
    opm/simulators/geochemistry/Core/GeoChemPhases.cpp
    opm/simulators/geochemistry/Core/GeoChemSolutionsManager.cpp
    opm/simulators/geochemistry/Core/MineralKinetics.cpp
    opm/simulators/geochemistry/Core/SimulationSchedule.cpp
    opm/simulators/geochemistry/Core/SurfaceChemistry.cpp
    opm/simulators/geochemistry/Database/ChemDatabaseHKF_aq.cpp
    opm/simulators/geochemistry/Database/ChemDatabaseHKF_basis.cpp
    opm/simulators/geochemistry/Database/ChemDatabaseHKF_buffer.cpp
    opm/simulators/geochemistry/Database/ChemDatabaseHKF_surface.cpp
    opm/simulators/geochemistry/IO/FileFunctions.cpp
    opm/simulators/geochemistry/IO/GeoLogger.cpp
    opm/simulators/geochemistry/IO/input_m.cpp
    opm/simulators/geochemistry/IO/InputFileHandler.cpp
    opm/simulators/geochemistry/IO/InputReader.cpp
    opm/simulators/geochemistry/IO/NetCDF.cpp
    opm/simulators/geochemistry/IO/ParseChemistry.cpp
    opm/simulators/geochemistry/IO/ParseUtil.cpp
    opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.cpp
    opm/simulators/geochemistry/SplayTree/SplayTree.cpp
    opm/simulators/geochemistry/SplayTree/SplayTreeKeyCreator.cpp
    opm/simulators/geochemistry/SplayTree/TreeNode.cpp
    opm/simulators/geochemistry/Thermo/eps_JN.cpp
    opm/simulators/geochemistry/Thermo/hkf.cpp
    opm/simulators/geochemistry/Thermo/ions.cpp
    opm/simulators/geochemistry/Thermo/ThermoTable.cpp
    opm/simulators/geochemistry/Thermo/thermodata.cpp
    opm/simulators/geochemistry/Thermo/water.cpp
    opm/simulators/geochemistry/Utility/HelpfulStringMethods.cpp
)

set(GEOCHEM_OPM_SOURCE_FILES
    opm/models/io/vtkgeochemistryparams.cpp
    opm/simulators/geochemistry/OpmGeoChemInterface.cpp
)

set(GEOCHEM_CORE_PUBLIC_HEADER_FILES
    opm/simulators/geochemistry/GeoChemIF.h
    opm/simulators/geochemistry/StandaloneSolvers.hpp
    opm/simulators/geochemistry/Common/ChemGlobal.h
    opm/simulators/geochemistry/Common/Constants.hpp
    opm/simulators/geochemistry/Common/CustomExceptions.hpp
    opm/simulators/geochemistry/Common/Enums.hpp
    opm/simulators/geochemistry/Common/SerializeForTesting.hpp
    opm/simulators/geochemistry/Core/ChemBasVec.h
    opm/simulators/geochemistry/Core/ChemGCSolver.h
    opm/simulators/geochemistry/Core/ChemInitChem.h
    opm/simulators/geochemistry/Core/ChemParam.h
    opm/simulators/geochemistry/Core/ChemState.hpp
    opm/simulators/geochemistry/Core/ChemTable.h
    opm/simulators/geochemistry/Core/GeoChemPhases.hpp
    opm/simulators/geochemistry/Core/GeoChemSolutionsManager.hpp
    opm/simulators/geochemistry/Core/MineralKinetics.hpp
    opm/simulators/geochemistry/Core/SimulationSchedule.hpp
    opm/simulators/geochemistry/Core/SurfaceChemistry.hpp
    opm/simulators/geochemistry/Database/ChemDatabaseHKF.hpp
    opm/simulators/geochemistry/Extra/DiffuseLayerIntegrand.hpp
    opm/simulators/geochemistry/IO/ErrorMsg.hpp
    opm/simulators/geochemistry/IO/fmt_include.hpp
    opm/simulators/geochemistry/IO/FileFunctions.hpp
    opm/simulators/geochemistry/IO/GeoLogger.hpp
    opm/simulators/geochemistry/IO/input_m.h
    opm/simulators/geochemistry/IO/InputFileHandler.hpp
    opm/simulators/geochemistry/IO/InputReader.hpp
    opm/simulators/geochemistry/IO/NetCDF.hpp
    opm/simulators/geochemistry/IO/ParseChemistry.hpp
    opm/simulators/geochemistry/IO/ParseUtil.hpp
    opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp
    opm/simulators/geochemistry/Numerical/BrentMethod.hpp
    opm/simulators/geochemistry/Numerical/GaussLegendre.hpp
    opm/simulators/geochemistry/Numerical/GaussLegendrePrecomputedTables.hpp
    opm/simulators/geochemistry/Numerical/NumericalErrorPolicy.hpp
    opm/simulators/geochemistry/Numerical/NumericalHelperFunctions.hpp
    opm/simulators/geochemistry/Numerical/RootSolvers.hpp
    opm/simulators/geochemistry/Numerical/SecantMethod.hpp
    opm/simulators/geochemistry/SplayTree/SplayTree.h
    opm/simulators/geochemistry/SplayTree/SplayTreeKeyCreator.hpp
    opm/simulators/geochemistry/SplayTree/TreeNode.h
    opm/simulators/geochemistry/Thermo/eps_JN.h
    opm/simulators/geochemistry/Thermo/hkf.h
    opm/simulators/geochemistry/Thermo/ions.h
    opm/simulators/geochemistry/Thermo/ThermoTable.h
    opm/simulators/geochemistry/Thermo/thermodata.h
    opm/simulators/geochemistry/Thermo/water.h
    opm/simulators/geochemistry/Utility/HelperMacros.hpp
    opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp
    opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp
)

set(GEOCHEM_OPM_PUBLIC_HEADER_FILES
    opm/models/io/vtkgeochemistrymodule.hpp
    opm/models/io/vtkgeochemistryparams.hpp
    opm/simulators/flow/FlowProblemGeochemistry.hpp
    opm/simulators/geochemistry/OpmGeoChemInterface.hpp
)

set(GEOCHEM_STANDALONE_EXAMPLE_SOURCE_FILES
    examples/GeoChemX.cpp
    examples/ThermoX.cpp
)

# OPM's build helpers still expect these legacy variables.
set(MAIN_SOURCE_FILES
    ${GEOCHEM_CORE_SOURCE_FILES}
    ${GEOCHEM_OPM_SOURCE_FILES}
)

set(PUBLIC_HEADER_FILES
    ${GEOCHEM_CORE_PUBLIC_HEADER_FILES}
    ${GEOCHEM_OPM_PUBLIC_HEADER_FILES}
)
