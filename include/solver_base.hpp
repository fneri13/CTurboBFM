#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"
#include "advection_base.hpp"
#include "advection_jst.hpp"
#include "advection_roe.hpp"
#include "boundary_base.hpp"
#include "boundary_inviscid_wall.hpp"
#include "boundary_inlet_2d.hpp"
#include "boundary_outlet.hpp"
#include "boundary_inlet_supersonic.hpp"
#include "boundary_outlet_supersonic.hpp"
#include "boundary_outlet_radial_equilibrium.hpp"
#include "boundary_fake.hpp"
#include "boundary_outlet_throttle.hpp"
#include "boundary_transparent.hpp"
#include <memory>
#include <vector>
#include <array>


class SolverBase {
    
public:

    SolverBase(Config &config, Mesh &mesh);
    
    virtual ~SolverBase() {}

    void readBoundaryConditions() ;

    virtual void initializeSolutionArrays() = 0;

    virtual void solve() = 0;

    const std::array<int, 3> getStepMask(FluxDirection direction) const;

    FloatType getHubStaticPressure() const { 
        return _hubStaticPressure; 
    }

    virtual void checkConvergence(bool &exitLoop, bool &skip) const = 0;

    virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) = 0;

    /** fetch indices for a 2D boundary slice of the 3D problem structure */
    void getBoundarySliceIndices(
        BoundaryIndices boundaryIdx, 
        size_t &iStart, 
        size_t &iLast, 
        size_t &jStart, 
        size_t &jLast, 
        size_t &kStart, 
        size_t &kLast) const;

protected:
    const Config& _config;
    Mesh& _mesh;
    size_t _nDimensions {0}, _nPointsI {0}, _nPointsJ {0}, _nPointsK {0};
    Topology _topology;
    
    Matrix3D<FloatType> _timeStep;
    std::vector<FloatType> _time;
    
    std::unique_ptr<FluidBase> _fluid;
    FluidModel _fluidModel = FluidModel::IDEAL;
    std::unique_ptr<AdvectionBase> _advection;
    
    std::map<BoundaryIndices, BoundaryType> _boundaryTypes;
    std::map<BoundaryIndices, std::vector<FloatType>> _boundaryValues;
    std::map<BoundaryIndices, Vector3D> _boundaryVelocities;
    std::map<BoundaryIndices, std::unique_ptr<BoundaryBase>> _boundaryConditions;
    FloatType _hubStaticPressure;
    std::vector<FloatType> _radialProfilePressure; 
    std::vector<FloatType> _radialProfileRadialCoords; 
    std::string _inlet2DfilePath{""}; 

    std::map<BoundaryIndices, FloatType> _massFlows;
    std::map<TurboPerformance, std::vector<FloatType>> _turboPerformance; 
    std::vector<std::map<MonitorOutputFields, std::vector<FloatType>>> _monitorPoints; 

    size_t _residualsDropConvergence {16};

};