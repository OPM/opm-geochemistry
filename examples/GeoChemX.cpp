/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <utility>

#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>
#include <opm/simulators/geochemistry/IO/InputFileHandler.hpp>
#include <opm/simulators/geochemistry/IO/GeoLogger.hpp>

int main(int argc, char* argv[])
{
    if (argc >= 2 && strcmp(argv[1], "EQSOLVER") == 0)
    {
        if (argc == 2)
        {
            std::cout << "Need to specify an input file for the equilibrium geochemical solver.\n";
        }
        else
        {
            auto inp = InputFileHandler(argv[2]);
            EquilibriumSolver equilibriumSolver;
            equilibriumSolver.solve(inp.getName(), inp.getStream());
        }
    }
    else if (argc >= 2 && strcmp(argv[1], "TRANSPORT") == 0)
    {
        if (argc == 2)
        {
            std::cout << "Need to specify an input file for the 1D transport solver.\n";
        }
        else
        {
            auto inp = InputFileHandler(argv[2]);
            OneDimensionalTransportSolver transportSolver;
            transportSolver.solve(inp.getName(), inp.getStream());
        }
    }
    else if (argc >= 2 && strcmp(argv[1], "THERMOTABLE") == 0)
    {
        if (argc == 2)
        {
            std::cout << "Need to specify an input file for the thermo table solver.\n";
        }
        else
        {
            auto inp = InputFileHandler(argv[2]);
            ThermoTableSolver thermoTableSolver;
            std::cout << thermoTableSolver.solve(inp.getName(), inp.getStream());
        }
    }
    else{
        std::cout << "No commands were given ... do not bother me unnecessarily\n!";
    }

    return 0;
}
