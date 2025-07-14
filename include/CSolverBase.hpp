#pragma once

#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"
#include "CFluid.hpp"
#include "CAdvectionSchemeBase.hpp"
#include "CAdvectionSchemeJST.hpp"
#include "CAdvectionSchemeROE.hpp"
#include "CBoundaryConditionBase.hpp"
#include "CBoundaryConditionEulerWall.hpp"
#include "CBoundaryConditionInlet.hpp"
#include "CBoundaryConditionOutlet.hpp"
#include "CBoundaryConditionInletSupersonic.hpp"
#include "CBoundaryConditionOutletSupersonic.hpp"
#include "CBoundaryConditionWedge.hpp"
#include "CBoundaryConditionRadialEquilibrium.hpp"
#include "CBoundaryConditionPeriodic.hpp"
#include "CBoundaryConditionThrottle.hpp"
#include <memory>


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

    protected:
        const Config& _config;
        
        const CMesh& _mesh;
        
        size_t _nDimensions {0}, _nPointsI, _nPointsJ, _nPointsK;
        
        Matrix3D<FloatType> _timeStep;

        std::vector<FloatType> _time;
        
        std::unique_ptr<CFluid> _fluid;

        std::unique_ptr<CAdvectionSchemeBase> _advectionScheme;
        
        std::map<BoundaryIndices, FloatType> _massFlows;
        
        std::map<BoundaryIndices, BoundaryType> _boundaryTypes;

        std::map<BoundaryIndices, std::vector<FloatType>> _boundaryValues;

        std::map<BoundaryIndices, std::unique_ptr<CBoundaryConditionBase>> _boundaryConditions;

        Topology _topology;

        FloatType _hubStaticPressure;

        std::vector<FloatType> _radialProfilePressure;
        
        std::vector<FloatType> _radialProfileCoords;

        std::map<TurboPerformance, std::vector<FloatType>> _turboPerformance; // map of turbo performance>

        std::vector<std::map<MonitorOutputFields, std::vector<FloatType>>> _monitorPoints; // vector of maps of monitor points data

        

};