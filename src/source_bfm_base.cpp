#include "source_bfm_base.hpp"


SourceBFMBase::SourceBFMBase(const Config &config, const FluidBase &fluid, const Mesh &mesh) 
        : _config(config), _fluid(fluid), _mesh(mesh) {

            _isBlockageActive = _config.isBlockageActive();
            _isPerturbationActive = _config.isPerturbationBodyForceActive();
            if (_isPerturbationActive){
                std::vector<FloatType> perturbationCenterCoords = _config.getPerturbationCenter();        
                _perturbationCenter = {
                    perturbationCenterCoords[0], 
                    perturbationCenterCoords[1], 
                    perturbationCenterCoords[2]};
                _perturbationRadius = _config.getPerturbationRadialExtension();
                _perturbationScaling = _config.getPerturbationScalingFactor();
                _perturbationTimeStart = _config.getPerturbationStartTime();
                _perturbationTimeDuration = _config.getPerturbationTimeDuration();
            }
  
        };


StateVector SourceBFMBase::computeTotalSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    Matrix3D<FloatType> &deviationAngle, 
    FloatType &timePhysical,
    FlowSolution &conservativeVars,
    Matrix3D<FloatType> &timestep) {

    StateVector blockageSource({0,0,0,0,0});
    if (_isBlockageActive){
        blockageSource = computeBlockageSource(i, j, k, primitive);
    }

    FloatType numberBlades = _mesh.getInputFields(InputField::NUMBER_BLADES, i, j, k);
    if (numberBlades == 0){
        return blockageSource;
    }

    FloatType dt = timestep(i, j, k);
    StateVector bodyForceSource = computeBodyForceSource(
        i, j, k, 
        primitive, 
        inviscidForce, 
        viscousForce, 
        conservativeVars,
        dt);
    
    FloatType perturbFactor = computePerturbationFactor(i, j, k, timePhysical);
    
    // update deviation angle field, needed for output
    deviationAngle(i, j, k) = _deviationAngle; 

    return blockageSource + bodyForceSource * (1.0 + perturbFactor);
}


StateVector SourceBFMBase::computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive) {

    Vector3D velocity({primitive[1], primitive[2], primitive[3]});
    FloatType totalEnthalpy = _fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType blockage = _mesh.getInputFields(InputField::BLOCKAGE, i, j, k);
    Vector3D blockageGrad = _mesh.getInputFieldsGradient(InputField::BLOCKAGE, i, j, k);
    FloatType volume = _mesh.getVolume(i, j, k);
    
    StateVector source({0,0,0,0,0});
    FloatType commonTerm = -1.0 / blockage * primitive[0] * velocity.dot(blockageGrad);
    source[0] = commonTerm;
    source[1] = commonTerm * velocity.x();
    source[2] = commonTerm * velocity.y();
    source[3] = commonTerm * velocity.z();
    source[4] = commonTerm * totalEnthalpy;

    return source*volume;
}


StateVector SourceBFMBase::computeBodyForceSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    FlowSolution &conservativeVars,
    FloatType &dt) {

    inviscidForce(i, j, k) = {0, 0, 0};
    viscousForce(i, j, k) = {0, 0, 0};
    return StateVector({0, 0, 0, 0, 0});
}


void SourceBFMBase::computeFlowState(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    FlowSolution &conservativeVars) {
    
    _radius = _mesh.getRadius(i, j, k);
    _theta = _mesh.getTheta(i, j, k);
    
    _velCartesian = {primitive[1], primitive[2], primitive[3]};
    _velCylindrical = computeCylindricalComponentsFromCartesian(_velCartesian, _theta);
    
    _omega = _mesh.getInputFields(InputField::RPM, i, j, k) * 2 * M_PI / 60;
    FloatType scalingOmega = _config.getRotationalSpeedScalingFactor();
    _omega *= scalingOmega;
    _dragVelCylindrical = {0, 0, _omega * _radius};
    _relVelCylindric = _velCylindrical - _dragVelCylindrical;
    _relVelCartesian = computeCartesianComponentsFromCylindrical(_relVelCylindric, _theta);
    _velMeridional = std::sqrt(_velCylindrical.x() * _velCylindrical.x() + 
                                _velCylindrical.y() * _velCylindrical.y());

    _bladesNumber = _mesh.getInputFields(InputField::NUMBER_BLADES, i, j, k);
    _tangentialPitch = 2.0 * M_PI * _radius / _bladesNumber;
    _normalCamberAxial = _mesh.getInputFields(InputField::NORMAL_AXIAL, i, j, k);
    _normalCamberRadial = _mesh.getInputFields(InputField::NORMAL_RADIAL, i, j, k);
    _normalCamberTangential = _mesh.getInputFields(InputField::NORMAL_TANGENTIAL, i, j, k);
    _normalCamberCylindric = {_normalCamberAxial, _normalCamberRadial, _normalCamberTangential};
    _blockage = _mesh.getInputFields(InputField::BLOCKAGE, i, j, k);
    _isBladePresent = _mesh.getInputFields(InputField::BLADE_PRESENT, i, j, k);
    
    _leanAngle = _mesh.getInputFields(InputField::BLADE_LEAN_ANGLE, i, j, k);
    _gasPathAngle = _mesh.getInputFields(InputField::BLADE_GAS_PATH_ANGLE, i, j, k);
    _metalAngle = _mesh.getInputFields(InputField::BLADE_METAL_ANGLE, i, j, k);
    _flowAngle = std::atan2(_relVelCylindric.z(), _velMeridional);
    _deviationAngle = computeDeviationAngle(_relVelCylindric, _normalCamberCylindric);
    
    // CONVENTION: deviation is positive when the blade should push the flow (flow is currently under-turned)
    if (_omega > 0.0){
        _deviationAngle = - _flowAngle + _metalAngle;
    }
    else {
        _deviationAngle = + _flowAngle - _metalAngle;
    }

    // CONVENTION: positive inviscid force means a positive force directed from suction to pressure side
    _inviscidForceDirCylindrical = computeInviscidForceDirection(_relVelCylindric, _normalCamberCylindric);
    _inviscidForceDirCartesian = computeCartesianComponentsFromCylindrical(_inviscidForceDirCylindrical, _theta);
    _viscousForceDirCylindrical = - _relVelCylindric.normalized();
    _viscousForceDirCartesian = - _relVelCartesian.normalized();

}


FloatType SourceBFMBase::computeDeviationAngle(Vector3D relativeVelocityCyl, Vector3D normalCamber){
    Vector3D normal = normalCamber / normalCamber.magnitude();
    Vector3D velMeridional = {relativeVelocityCyl.x(), relativeVelocityCyl.y(), 0.0};
    Vector3D spanwiseDir = {-velMeridional.y(), velMeridional.x(), 0.0};
    spanwiseDir = spanwiseDir / spanwiseDir.magnitude();

    // blade-to-blade velocity obtained removing the spanwise component
    Vector3D velInPlane = relativeVelocityCyl - spanwiseDir * relativeVelocityCyl.dot(spanwiseDir);
    FloatType normalVel = velInPlane.dot(normal);

    // the deviation angle is positive when the blade should push, and negative when the blade should pull
    FloatType deviationAngle = -std::asin(normalVel / velInPlane.magnitude());
    
    return deviationAngle;
}

Vector3D SourceBFMBase::computeInviscidForceDirection(const Vector3D& relativeVelocity, const Vector3D& normalCamber){

    // remove component along relative velocity
    Vector3D relVelDir = relativeVelocity.normalized();
    Vector3D normalDir = normalCamber.normalized();
    Vector3D versor = (normalDir - relVelDir * normalDir.dot(relVelDir)).normalized();

    // remove component out of blade-to-blade plane
    Vector3D velMeridional{relativeVelocity.x(), relativeVelocity.y(), 0.0};
    Vector3D spanDirection = {velMeridional.y(), -velMeridional.x(), 0.0};
    spanDirection = spanDirection.normalized();
    versor = (versor - spanDirection * versor.dot(spanDirection)).normalized();

    return versor;

}


FloatType SourceBFMBase::computeTangentialComponent(
    FloatType fAxial, 
    const Vector3D& relativeVelocityDirection, 
    const Vector3D& normalCamber){

    FloatType fTangential = (
        -relativeVelocityDirection.x() * fAxial -
         relativeVelocityDirection.y() * normalCamber.y()) /
        relativeVelocityDirection.z();
    return fTangential;
}


FloatType SourceBFMBase::computePerturbationFactor(size_t i, size_t j, size_t k, FloatType timePhysical) const{

    if (!_config.isPerturbationBodyForceActive()){
        return 0.0;
    }

    if (timePhysical < _perturbationTimeStart || timePhysical > _perturbationTimeStart + _perturbationTimeDuration){
        return 0.0;
    }
    
    FloatType variationMax = _perturbationScaling; 

    Vector3D point = _mesh.getVertex(i,j,k);
    FloatType distance = (point - _perturbationCenter).magnitude();

    if (distance > 3.0*_perturbationRadius){
        return 0.0;
    }

    FloatType spaceModulation = std::exp(-(distance*distance)/(_perturbationRadius*_perturbationRadius));

    FloatType peakInTime = _perturbationTimeStart + _perturbationTimeDuration/2.0;
    FloatType sigmaTime = _perturbationTimeDuration/3.0;
    FloatType timeModulation = std::exp( - (std::pow(timePhysical-peakInTime, 2)) / std::pow(sigmaTime, 2)); 

    FloatType perturbationVariation = variationMax * spaceModulation * timeModulation; 

    return perturbationVariation;
}
