#include "solver_base.hpp"
#include <iostream>

SolverBase::SolverBase(Config& config, Mesh& mesh)
    : _config(config), _mesh(mesh)
{
    _nDimensions = _mesh.getNumberDimensions();
    _nPointsI = _mesh.getNumberPointsI();
    _nPointsJ = _mesh.getNumberPointsJ();
    _nPointsK = _mesh.getNumberPointsK();

    _timeStep.resize(_nPointsI, _nPointsJ, _nPointsK);
    _time.push_back(0.0);

    _topology = _config.getTopology();

    _residualsDropConvergence = _config.getResidualsDropConvergence();  

    _fluidModel = _config.getFluidModel();
    if (_fluidModel == FluidModel::IDEAL){
        _fluid = std::make_unique<FluidIdeal>(_config.getFluidGamma(), _config.getFluidGasConstant());
    }
    else if (_fluidModel == FluidModel::REAL){
        throw std::runtime_error("Real fluid model not implemented yet.");
        // _fluid = std::make_unique<CFluidReal>(_config.getFluidRealDataFile());
    }
    else{
        throw std::runtime_error("Unsupported fluid model selected.");
    }

    ConvectionScheme advScheme = _config.getConvectionScheme();
    switch (advScheme)
    {
    case ConvectionScheme::JST:
        _advectionScheme = std::make_unique<AdvectionJst>(_config, *_fluid);
        break;
    case ConvectionScheme::ROE:
        _advectionScheme = std::make_unique<AdvectionRoe>(_config, *_fluid);
        break;
    default:
        throw std::runtime_error("Unsupported convection scheme selected.");
    }

    readBoundaryConditions();
    
}

void SolverBase::readBoundaryConditions(){
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
        else if (_boundaryTypes[bound] == BoundaryType::INLET_2D){
            _boundaryValues[bound] = std::vector<FloatType> {};
            _inlet2DfilePath = _config.getInlet2DfilePath();
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET || _boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC || _boundaryTypes[bound] == BoundaryType::THROTTLE){
            _boundaryValues[bound] = _config.getOutletBCValues();
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryValues[bound] = _config.getOutletBCValues();
            _hubStaticPressure = _boundaryValues[bound][0];
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _boundaryValues[bound] = _config.getPeriodicityInfo();
            if (!periodicityChecked) {
                FloatType angle = _config.getPeriodicityAngleRad();
                FloatType translation = _config.getPeriodicityTranslation();
                _mesh.checkPeriodicity(angle, translation);
                _mesh.setPeriodicMesh(angle, translation);
                periodicityChecked = true;
            }
        }
        else if (_boundaryTypes[bound] == BoundaryType::NO_SLIP_WALL){
            _boundaryVelocities[bound] = _config.getNoSlipWallVelocity(bound);
        }
        else {
            _boundaryValues[bound] = std::vector<FloatType> {};
        }
    }

    // instantiate boundary condition objects
    for (auto& bound : bounds) {
        if (_boundaryTypes[bound] == BoundaryType::WALL || _boundaryTypes[bound] == BoundaryType::NO_SLIP_WALL){
            _boundaryConditions[bound] = std::make_unique<BoundaryInviscidWall>(_config, _mesh, *_fluid, bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET){
            _boundaryConditions[bound] = std::make_unique<BoundaryInlet>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_2D){
            _boundaryConditions[bound] = std::make_unique<BoundaryInlet2D>(_config, _mesh, *_fluid, bound, _inlet2DfilePath);
        }
        else if (_boundaryTypes[bound] == BoundaryType::INLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryInletSupersonic>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutlet>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::RADIAL_EQUILIBRIUM){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletRadialEquilibrium>(_config, _mesh, *_fluid, bound, _radialProfilePressure);
        }
        else if (_boundaryTypes[bound] == BoundaryType::OUTLET_SUPERSONIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletSupersonic>(_config, _mesh, *_fluid, bound, _boundaryValues[bound]);
        }
        else if (_boundaryTypes[bound] == BoundaryType::THROTTLE){
            _boundaryConditions[bound] = std::make_unique<BoundaryOutletThrottle>(_config, _mesh, *_fluid, bound, _radialProfilePressure);
        }
        else if (_boundaryTypes[bound] == BoundaryType::WEDGE){
            _boundaryConditions[bound] = std::make_unique<BoundaryFake>(_config, _mesh, *_fluid, bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::PERIODIC){
            _boundaryConditions[bound] = std::make_unique<BoundaryFake>(_config, _mesh, *_fluid, bound);
        }
        else if (_boundaryTypes[bound] == BoundaryType::TRANSPARENT){
            assert (_mesh.getNumberPointsJ() == 1 && "Transparent boundary only valid for 1D meshes.");
            assert (_mesh.getNumberPointsK() == 1 && "Transparent boundary only valid for 1D meshes.");
            _boundaryConditions[bound] = std::make_unique<BoundaryTransparent>(_config, _mesh, *_fluid, bound, *_advectionScheme);
        }
    }
}


const std::array<int, 3> SolverBase::getStepMask(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return {1, 0, 0};
        case FluxDirection::J: return {0, 1, 0};
        case FluxDirection::K: return {0, 0, 1};
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}


void SolverBase::getBoundarySliceIndices(BoundaryIndices boundaryIdx, size_t &iStart, size_t &iLast, size_t &jStart, size_t &jLast, size_t &kStart, size_t &kLast) const{
    iStart=0, 
    iLast=_nPointsI, 
    jStart=0, 
    jLast=_nPointsJ, 
    kStart=0, 
    kLast=_nPointsK;
    
    switch (boundaryIdx)
    {
    case BoundaryIndices::I_START:
        iLast = 1;
        break;
    case BoundaryIndices::I_END:
        iStart = _nPointsI-1;
        break;
    case BoundaryIndices::J_START:
        jLast = 1;
        break;
    case BoundaryIndices::J_END:
        jStart = _nPointsJ-1;
        break;
    case BoundaryIndices::K_START:
        kLast = 1;
        break;
    case BoundaryIndices::K_END:
        kStart = _nPointsK-1;
        break;
    }
}