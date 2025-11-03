#pragma once

#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"
#include "CFluidBase.hpp"
#include "CFluidIdeal.hpp"
#include "CAdvectionSchemeBase.hpp"
#include "CAdvectionSchemeJST.hpp"
#include "CAdvectionSchemeROE.hpp"
#include "CBoundaryConditionBase.hpp"
#include "CBoundaryConditionEulerWall.hpp"
#include "CBoundaryConditionInlet.hpp"
#include "CBoundaryConditionInlet2D.hpp"
#include "CBoundaryConditionOutlet.hpp"
#include "CBoundaryConditionInletSupersonic.hpp"
#include "CBoundaryConditionOutletSupersonic.hpp"
#include "CBoundaryConditionWedge.hpp"
#include "CBoundaryConditionRadialEquilibrium.hpp"
#include "CBoundaryConditionPeriodic.hpp"
#include "CBoundaryConditionThrottle.hpp"
#include <memory>
#include <vector>
#include <array>


/** 
  *  \brief     Class handling base solver capabilities.
  *  \details   It computes the finite volume mesh.
  *  \author    Francesco Neri
  */
class CSolverBase {
    
    public:

        CSolverBase(Config &config, CMesh &mesh);
        
        virtual ~CSolverBase() {}

        void readBoundaryConditions();

        virtual void initializeSolutionArrays() = 0;

        virtual void solve() = 0;

        const std::array<int, 3> getStepMask(FluxDirection direction) const;

        FloatType getHubStaticPressure() const { return _hubStaticPressure; }

        virtual void checkConvergence(bool &exitLoop, bool &skip) const = 0;

        virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) = 0;

    protected:
        const Config& _config;
        
        CMesh& _mesh;
        
        size_t _nDimensions {0}, _nPointsI, _nPointsJ, _nPointsK;
        
        Matrix3D<FloatType> _timeStep;

        std::vector<FloatType> _time;
        
        std::unique_ptr<CFluidBase> _fluid;

        std::unique_ptr<CAdvectionSchemeBase> _advectionScheme;
        
        std::map<BoundaryIndices, FloatType> _massFlows;
        
        std::map<BoundaryIndices, BoundaryType> _boundaryTypes;

        std::map<BoundaryIndices, std::vector<FloatType>> _boundaryValues;

        std::map<BoundaryIndices, Vector3D> _boundaryVelocities;

        std::map<BoundaryIndices, std::unique_ptr<CBoundaryConditionBase>> _boundaryConditions;

        Topology _topology;

        FloatType _hubStaticPressure;

        std::vector<FloatType> _radialProfilePressure; // to be used from BC with radial equilibrium
        
        std::vector<FloatType> _radialProfileCoords; // to be used from BC with radial equilibrium

        std::map<TurboPerformance, std::vector<FloatType>> _turboPerformance; // map of turbo performance>

        std::vector<std::map<MonitorOutputFields, std::vector<FloatType>>> _monitorPoints; // vector of maps of monitor points data

        size_t _residualsDropConvergence = 16;

        FluidModel _fluidModel = FluidModel::IDEAL;

        std::string _inlet2DfilePath{""}; // name of the 2D inlet file

        /* fetch indices for a 2D boundary slice depdending on boundary index */
        void fetchBoundarySliceIndices(BoundaryIndices boundaryIdx, size_t &iStart, size_t &iLast, size_t &jStart, size_t &jLast, size_t &kStart, size_t &kLast) const;

        

};