#pragma once
#include "boundary_base.hpp"

class BoundaryInviscidWall : public BoundaryBase {
    
public:
    
    /**
     * @brief Constructs the object with given references.
     * @param config The configuration object.
     * @param mesh The mesh object.
     * @param fluid The fluid object.
     * @param boundIndex The boundary index.
     */
    BoundaryInviscidWall(
        const Config &config, 
        const Mesh &mesh, 
        const FluidBase &fluid, 
        BoundaryIndices boundIndex)
        : BoundaryBase(config, mesh, fluid, boundIndex) {};
        
    virtual ~BoundaryInviscidWall() = default;

    virtual StateVector computeBoundaryFlux(
        const StateVector& internalConservative, 
        const Vector3D& surface, 
        const Vector3D& midPoint, 
        const std::array<size_t, 3>& indices, 
        const FlowSolution& flowSolution, 
        const size_t& iterCounter) override;
              
};