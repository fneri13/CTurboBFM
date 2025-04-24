#include "../include/CSolverBase.hpp"
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
        _advectionScheme = std::make_unique<CJSTScheme>(*_fluid);
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

    // read the boundaries values
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::INLET){
            _boundaryValues[bound] = _config.getInletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET){
            _boundaryValues[bound] = _config.getOutletBCValues();
        }
        else {
            _boundaryValues[bound] = std::vector<FloatType> {};
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

