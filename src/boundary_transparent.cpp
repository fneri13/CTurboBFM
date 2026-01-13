#include "boundary_transparent.hpp"
#include "math_utils.hpp"

BoundaryTransparent::BoundaryTransparent(const Config &config, const Mesh &mesh, FluidBase &fluid, BoundaryIndices boundIndex, AdvectionBase &advScheme)
    : BoundaryBase(config, mesh, fluid, boundIndex), _advScheme(advScheme) {
    }


    StateVector BoundaryTransparent::computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) {
        StateVector Ul = internalConservative;
        StateVector Ur = internalConservative;
        StateVector flux = _advScheme.computeFlux(Ul, Ul, Ur, Ur, surface);
        return flux;
    }