#include "CSolverBase.hpp"
#include <iostream>

CSolverBase::CSolverBase(Config& config, CMesh& mesh)
    : _config(config), _mesh(mesh)
{
    _nDimensions = _mesh.getNumberDimensions();
    _nPointsI = _mesh.getNumberPointsI();
    _nPointsJ = _mesh.getNumberPointsJ();
    _nPointsK = _mesh.getNumberPointsK();
    _timeStep.resize(_nPointsI, _nPointsJ, _nPointsK);

    _topology = _config.getTopology();

    _fluid = std::make_unique<CFluid>(_config.getFluidGamma(), _config.getFluidGasConstant());

    ConvectionScheme advScheme = _config.getConvectionScheme();

    switch (advScheme)
    {
    case ConvectionScheme::JST:
        _advectionScheme = std::make_unique<CAdvectionSchemeJST>(*_fluid);
        break;

    default:
        throw std::runtime_error("Unsupported convection scheme selected.");
    }

    readBoundaryConditions();
    
}

void CSolverBase::readBoundaryConditions(){
    std::array<BoundaryIndices, 6> bounds = {BoundaryIndices::I_START,
        BoundaryIndices::I_END,
        BoundaryIndices::J_START,
        BoundaryIndices::J_END,
        BoundaryIndices::K_START,
        BoundaryIndices::K_END};
    
    // read the boundaries type
    for (auto& bound : bounds) {
        _boundaryTypes[bound] = _config.getBoundaryType(bound);
    }

    bool periodicityChecked = false;

    // read the boundaries values
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::INLET || _boundaryTypes[bound] == BoundaryType::INLET_SUPERSONIC) {
            _boundaryValues[bound] = _config.getInletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET || _boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC){
            _boundaryValues[bound] = _config.getOutletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryValues[bound] = _config.getOutletBCValues();
            _hubStaticPressure = _boundaryValues[bound][0];
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _boundaryValues[bound].push_back(_config.getPeriodicityAngleRad());
            if (!periodicityChecked) {
                _mesh.checkPeriodicity(_boundaryValues[bound].back());
                periodicityChecked = true;
            }
        }
        else {
            _boundaryValues[bound] = std::vector<FloatType> {};
        }
    }

    // instantiate boundary condition objects
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::WALL){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionEulerWall>(_config, _mesh, *_fluid, bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionInlet>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionInletSupersonic>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionOutlet>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionRadialEquilibrium>(_config, _mesh, *_fluid, bound, _radialProfilePressure);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionOutletSupersonic>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::WEDGE){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionWedge>(_config, _mesh, *_fluid, bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _boundaryConditions[bound] = std::make_unique<CBoundaryConditionPeriodic>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
    }
}


const std::array<int, 3> CSolverBase::getStepMask(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return {1, 0, 0};
        case FluxDirection::J: return {0, 1, 0};
        case FluxDirection::K: return {0, 0, 1};
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}

