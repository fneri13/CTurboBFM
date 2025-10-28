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
    // properties of internal point
    StateVector primitive = getEulerPrimitiveFromConservative(internalConservative);
    Vector3D velocityInt({primitive[1], primitive[2], primitive[3]});
    FloatType soundSpeedInt = _fluid.computeSoundSpeed_rho_u_et(primitive[0], velocityInt, primitive[4]);
    FloatType totEnthalpyInt = _fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocityInt, primitive[4]);
    FloatType Jm = - velocityInt.magnitude() + 2*soundSpeedInt / (_fluid.getGamma() - 1);

    // Solve the quadratic equation for the speed of sound
    FloatType alpha = 1.0 / (_fluid.getGamma() - 1.0) + 2.0 / std::pow((_fluid.getGamma() - 1.0), 2);
    FloatType beta = -2.0 * Jm / (_fluid.getGamma() - 1.0);
    FloatType zeta = 0.5 * Jm * Jm - totEnthalpyInt;
    FloatType soundSpeedBound = std::max((-beta + std::sqrt(beta*beta - 4.0*alpha*zeta))/2.0/alpha,
                                         (-beta - std::sqrt(beta*beta - 4.0*alpha*zeta))/2.0/alpha);



    // reconstruct the boundary state, taking info from the boundary                                      
    FloatType velocityBoundMag = 2.0*soundSpeedBound / (_fluid.getGamma() - 1.0) - Jm;
    
    // flow direction from input file
    FloatType nx = _inputGrid.getField(FieldNames::INLET_NX, indices[0], indices[1], indices[2]);
    FloatType ny = _inputGrid.getField(FieldNames::INLET_NY, indices[0], indices[1], indices[2]);
    FloatType nz = _inputGrid.getField(FieldNames::INLET_NZ, indices[0], indices[1], indices[2]);
    Vector3D flowDirection(nx, ny, nz);
    if (_referenceFrame == ReferenceFrame::CYLINDRICAL){
        FloatType theta = _mesh.getTheta(indices[0], indices[1], indices[2]);
        flowDirection = computeCartesianVectorFromCylindrical(flowDirection, theta);
    }
    flowDirection /= flowDirection.magnitude();

    // boundary values from input file
    FloatType totalPressure = _inputGrid.getField(FieldNames::TOTAL_PRESSURE, indices[0], indices[1], indices[2]);
    FloatType totalTemperature = _inputGrid.getField(FieldNames::TOTAL_TEMPERATURE, indices[0], indices[1], indices[2]);
    
    FloatType normalMachBound = velocityBoundMag / soundSpeedBound;
    FloatType pressureBound = _fluid.computeStaticPressure_pt_M(totalPressure, normalMachBound);
    FloatType temperatureBound = _fluid.computeStaticTemperature_Tt_M(totalTemperature, normalMachBound);
    FloatType densityBound = _fluid.computeDensity_p_T(pressureBound, temperatureBound);
    FloatType energyBound = _fluid.computeStaticEnergy_p_rho(pressureBound, densityBound);
    Vector3D velocityBound = flowDirection * velocityBoundMag;
    FloatType totEnergyBound = energyBound + 0.5 * velocityBound.dot(velocityBound);

    // compute boundary flux
    StateVector primitiveBoundary({densityBound, velocityBound.x(), velocityBound.y(), velocityBound.z(), totEnergyBound});
    StateVector flux = computeEulerFluxFromPrimitive(primitiveBoundary, surface, _fluid);
    return flux;
}