#include "../../include/Config.hpp"
#include <iostream>

using namespace std;

int main() {
    Config config("input.ini");

    cout << "Grid file path: " << config.gridFilePath() << "\n";
    cout << "Restart solution: " << (config.restartSolution() ? "true" : "false") << "\n";
    cout << "Restart solution file path: " << config.restartSolutionFilepath() << "\n";
    cout << "Is BFM active: " << (config.isBFMActive() ? "true" : "false") << "\n";
    cout << "Is blockage active: " << (config.isBlockageActive() ? "true" : "false") << "\n";
    cout << "BFM model: " << static_cast<int>(config.getBFMModel()) << "\n";
    cout << "Kind solver: " << static_cast<int>(config.getKindSolver()) << "\n";
    cout << "Topology: " << static_cast<int>(config.getTopology()) << "\n";
    cout << "Fluid gamma: " << config.getFluidGamma() << "\n";
    cout << "Fluid gas constant: " << config.getFluidGasConstant() << "\n";
    cout << "Fluid kinematic viscosity: " << config.getFluidKinematicViscosity() << "\n";
    
    std::vector<BoundaryType> boundaryTypeI = config.getBoundaryTypeI();
    for (const auto& type : boundaryTypeI) {
        std::cout << "Boundary type I: " << static_cast<int>(type) << "\n";
    }

    std::vector<BoundaryType> boundaryTypeJ = config.getBoundaryTypeJ();
    for (const auto& type : boundaryTypeJ) {
        std::cout << "Boundary type J: " << static_cast<int>(type) << "\n";
    }

    std::vector<BoundaryType> boundaryTypeK = config.getBoundaryTypeK();
    for (const auto& type : boundaryTypeK) {
        std::cout << "Boundary type K: " << static_cast<int>(type) << "\n";
    }

    auto inletBCValues = config.getInletBCValues();
    cout << "Inlet BC values: ";
    for (const auto& value : inletBCValues) {
        cout << value << " ";
    }
    cout << "\n";

    auto outletBCValues = config.getOutletBCValues();
    cout << "Outlet BC values: ";
    for (const auto& value : outletBCValues) {
        cout << value << " ";
    }
    cout << "\n";

    cout << "Initial Mach number: " << config.getInitMachNumber() << "\n";
    cout << "Initial pressure: " << config.getInitPressure() << "\n";
    cout << "Initial temperature: " << config.getInitTemperature() << "\n";
    auto initDirection = config.getInitDirection();
    cout << "Initial direction: ";
    for (const auto& value : initDirection) {
        cout << value << " ";
    }
    cout << "\n";
    cout << "CFL: " << config.getCFL() << "\n";
    cout << "Max iterations: " << config.getMaxIterations() << "\n";
    cout << "Save unsteady solution: " << (config.saveUnsteadySolution() ? "true" : "false") << "\n";
    cout << "Time integration: " << static_cast<int>(config.getTimeIntegration()) << "\n";
    cout << "Time step method: " << static_cast<int>(config.getTimeStepMethod()) << "\n";
    cout << "Advection scheme: " << static_cast<int>(config.getConvectionScheme()) << "\n";

    return 0;
}