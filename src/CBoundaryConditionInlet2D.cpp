#include "CBoundaryConditionInlet2D.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionInlet2D::CBoundaryConditionInlet2D(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex, std::string inletFile)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        
        _inletFilePath = inletFile;
        
        _referenceFrame = config.getInletReferenceFrame();
        
        // read input file and check that is compatible with the mesh
        _inputGrid = CInput(_inletFilePath);
        
        _inputNi = _inputGrid.getNumberPointsI();
        _inputNj = _inputGrid.getNumberPointsJ();
        _inputNk = _inputGrid.getNumberPointsK();
        assert(_inputNi == 1 && "The inlet 2D file must have only one node along i.");
        assert((_inputNj == _mesh.getNumberPointsJ() && _inputNk == _mesh.getNumberPointsK()) &&
            "The inlet 2D file and the mesh must share the same number of nodes along j and k.");

    }


StateVector CBoundaryConditionInlet2D::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
    
    // flow direction from input file
    FloatType nx = _inputGrid.getField(FieldNames::INLET_NX, indices[0], indices[1], indices[2]);
    FloatType ny = _inputGrid.getField(FieldNames::INLET_NY, indices[0], indices[1], indices[2]);
    FloatType nz = _inputGrid.getField(FieldNames::INLET_NZ, indices[0], indices[1], indices[2]);
    Vector3D flowDirection(nx, ny, nz);
    flowDirection /= flowDirection.magnitude();

    // boundary values from input file
    FloatType totalPressure = _inputGrid.getField(FieldNames::TOTAL_PRESSURE, indices[0], indices[1], indices[2]);
    FloatType totalTemperature = _inputGrid.getField(FieldNames::TOTAL_TEMPERATURE, indices[0], indices[1], indices[2]);
    
    StateVector flux = computeSubsonicsInletFlux(internalConservative, surface, totalPressure, totalTemperature, flowDirection);
    
    return flux;
}