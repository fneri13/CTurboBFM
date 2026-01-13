#include "boundary_outlet_supersonic.hpp"
#include "math_utils.hpp"

BoundaryOutletSupersonic::BoundaryOutletSupersonic(const Config &config, const Mesh &mesh, FluidBase &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues)
    : BoundaryBase(config, mesh, fluid, boundIndex) {
        _boundaryValues = inletValues;
    }


StateVector BoundaryOutletSupersonic::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
    auto primitive = getEulerPrimitiveFromConservative(internalConservative);
    auto flux = computeEulerFluxFromPrimitive(primitive, surface, _fluid);
    return flux;
}