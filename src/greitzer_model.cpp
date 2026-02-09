#include "greitzer_model.hpp"

GreitzerModel::GreitzerModel(
    const Config &config, 
    const FluidBase &fluid)
    : _fluid(fluid), _config(config)
{
    _plenumVolume = _config.getGreitzerPlenumVolume();
    _throttleCoefficient = _config.getGreitzerThrottleCoefficient();
    _deltaTime = _config.getFixedTimeStep();
    _fluidGamma = _fluid.getGamma();
    _fluidRConstant = _fluid.getRconstant();
}

void GreitzerModel::initializeState(
    FloatType plenumPressure, 
    FloatType plenumInletMassflow, 
    FloatType plenumOutletMassflow) {
    
    _plenumPressure.push_back(plenumPressure);
    _plenumInletMassflow.push_back(plenumInletMassflow);
    _plenumOutletMassflow.push_back(plenumOutletMassflow);
    _time.push_back(0.0);
}



FloatType GreitzerModel::computePlenumPressure(FloatType massFlow) {
    _plenumInletMassflow.push_back(massFlow);
    _time.push_back(_time.back() + _deltaTime);
    
    FloatType aPlenum = std::sqrt(_fluidGamma * _fluidRConstant * 350); // for now keep it simple, it will change a bit the dynamics

    FloatType newP = _plenumPressure.back() + _deltaTime*aPlenum*aPlenum / (_fluidGamma*_plenumVolume) * (
        _plenumInletMassflow.back() - _plenumOutletMassflow.back());
    _plenumPressure.push_back(newP);

    FloatType newM = std::sqrt((_plenumPressure.back()-101325) / _throttleCoefficient);
    _plenumOutletMassflow.push_back(newM);
    
    return _plenumPressure.back();
}