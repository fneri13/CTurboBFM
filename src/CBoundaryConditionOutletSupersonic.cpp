#include "CBoundaryConditionOutletSupersonic.hpp"
#include "commonFunctions.hpp"

CBoundaryConditionOutletSupersonic::CBoundaryConditionOutletSupersonic(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : CBoundaryConditionBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
    }


StateVector CBoundaryConditionOutletSupersonic::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution) {
    auto primitive = getEulerPrimitiveFromConservative(internalConservative);
    auto flux = computeEulerFluxFromPrimitive(primitive, surface, _fluid);
    return flux;
}