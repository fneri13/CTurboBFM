#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"


class BoundaryBase {
    
public:
    BoundaryBase(const Config &config, const Mesh &mesh, const FluidBase &fluid, BoundaryIndices boundIndex) 
        : _config(config), _mesh(mesh), _fluid(fluid), _boundaryIndex(boundIndex) {};
        
    virtual ~BoundaryBase() = default;

    /**
     * @brief Compute the boundary flux.
     * @param internalConservative The internal point conservative variables.
     * @param surface The surface normal vector.
     * @param midPoint The midpoint of the boundary face.
     * @param indices The indices (i,j,k) of the boundary face.
     * @param flowSolution The entire flow solution.
     */
    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) = 0;

protected:

    /** Common interface for subsonic inlet types. */
    StateVector computeSubsonicInletFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const FloatType& totPressure, 
        const FloatType& totTemperature, 
        const Vector3D& flowDirection);
    
protected:

    const Config& _config;
    const Mesh& _mesh;
    const FluidBase& _fluid;
    const BoundaryIndices _boundaryIndex;
        
};