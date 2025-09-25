#include "CSourceBFMBase.hpp"


CSourceBFMBase::CSourceBFMBase(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : _config(config), _fluid(fluid), _mesh(mesh) {

            // read data associated to body force perturbation
            _perturbationBodyForceActive = _config.isPerturbationBodyForceActive();
            if (_perturbationBodyForceActive){
                
                std::vector<FloatType> perturbationCenterCoords = _config.getPerturbationCenter();
                _perturbationCenter = {perturbationCenterCoords[0], perturbationCenterCoords[1], perturbationCenterCoords[2]};
                _perturbationRadius = _config.getPerturbationRadialExtension();
                _perturbationScalingFactor = _config.getPerturbationScalingFactor();
                _perturbationTimeStart = _config.getPerturbationStartTime();
                _perturbationTimeDuration = _config.getPerturbationTimeDuration();
                
            }
  
        };


StateVector CSourceBFMBase::computeSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, Matrix3D<FloatType> &deviationAngle, FloatType &timePhysical) {

    StateVector blockageSource({0,0,0,0,0});
    if (_config.isBlockageActive()){
        blockageSource = computeBlockageSource(i, j, k, primitive);
    }

    FloatType numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    if (numberBlades == 0){ // this is the layer upstream of leading edge, no blade present, only blockage contribution
        return blockageSource;
    }

    StateVector bodyForceSource = computeBodyForceSource(i, j, k, primitive, inviscidForce, viscousForce);
    FloatType bodyForceVariation = computeBodyForcePerturbationScaling(i, j, k, timePhysical);

    deviationAngle(i, j, k) = _deviationAngle; // update the deviation angle that has been computed from every bfm model, getting it ready for the output

    return blockageSource + bodyForceSource * (1.0 + bodyForceVariation);
}

StateVector CSourceBFMBase::computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive) {
    Vector3D velocity({primitive[1], primitive[2], primitive[3]});
    FloatType totalEnthalpy = _fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType blockage = _mesh.getInputFields(FieldNames::BLOCKAGE, i, j, k);
    Vector3D blockageGrad = _mesh.getInputFieldsGradient(FieldNames::BLOCKAGE, i, j, k);
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


StateVector CSourceBFMBase::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    inviscidForce(i, j, k) = {0, 0, 0};
    viscousForce(i, j, k) = {0, 0, 0};
    return StateVector({0, 0, 0, 0, 0});
}

void CSourceBFMBase::computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive){
    // geometrical considerations
    _point = _mesh.getVertex(i, j, k);
    _radius = std::sqrt(_point.z() * _point.z() + _point.y() * _point.y());
    _theta = std::atan2(_point.z(), _point.y());
    
    // velocity vector in both reference frame
    _velocityCartesian = {primitive[1], primitive[2], primitive[3]};
    _velocityCylindrical = computeCylindricalVectorFromCartesian(_velocityCartesian, _theta);
    
    // relative velocity vector in both reference frame
    _omega = _mesh.getInputFields(FieldNames::RPM, i, j, k) * 2 * M_PI / 60;
    _dragVelocityCylindrical = {0, 0, _omega * _radius};
    _relativeVelocityCylindric = _velocityCylindrical - _dragVelocityCylindrical;
    _relativeVelocityCartesian = computeCartesianVectorFromCylindrical(_relativeVelocityCylindric, _theta);
    _velMeridional = std::sqrt(_velocityCylindrical.x() * _velocityCylindrical.x() + 
                                _velocityCylindrical.y() * _velocityCylindrical.y());

    // blade properties
    _numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    _pitch = 2.0 * M_PI * _radius / _numberBlades;
    _normalCamberAxial = _mesh.getInputFields(FieldNames::NORMAL_AXIAL, i, j, k);
    _normalCamberRadial = _mesh.getInputFields(FieldNames::NORMAL_RADIAL, i, j, k);
    _normalCamberTangential = _mesh.getInputFields(FieldNames::NORMAL_TANGENTIAL, i, j, k);
    _normalCamberCylindric = {_normalCamberAxial, _normalCamberRadial, _normalCamberTangential};
    _blockage = _mesh.getInputFields(FieldNames::BLOCKAGE, i, j, k);
    _bladeIsPresent = _mesh.getInputFields(FieldNames::BLADE_PRESENT, i, j, k);
    
    // flow directions 
    
    _leanAngle = _mesh.getInputFields(FieldNames::BLADE_LEAN_ANGLE, i, j, k);
    _gaspathAngle = _mesh.getInputFields(FieldNames::BLADE_GAS_PATH_ANGLE, i, j, k);
    _metalAngle = _mesh.getInputFields(FieldNames::BLADE_METAL_ANGLE, i, j, k);
    _flowAngle = std::atan2(_relativeVelocityCylindric.z(), _velMeridional);
    _deviationAngle = computeDeviationAngle(_relativeVelocityCylindric, _normalCamberCylindric);
    // _deviationAngle = _flowAngle - _metalAngle;

    // compute inviscid force directions with the two methods
    _inviscidForceDirectionCylindrical = computeInviscidForceDirection(_relativeVelocityCylindric, _normalCamberCylindric);
    // _inviscidForceDirectionCylindrical = computeInviscidForceDirection(_relativeVelocityCylindric, _gaspathAngle, _flowAngle, _leanAngle);

    _inviscidForceDirectionCartesian = computeCartesianVectorFromCylindrical(_inviscidForceDirectionCylindrical, _theta);
    _viscousForceDirectionCylindrical = - _relativeVelocityCylindric.normalized();
    _viscousForceDirectionCartesian = - _relativeVelocityCartesian.normalized();
}


FloatType CSourceBFMBase::computeDeviationAngle(Vector3D relativeVelocity, Vector3D normalCamber){
    Vector3D normal = normalCamber / normalCamber.magnitude();
    FloatType normalVelocity = relativeVelocity.dot(normal);

    // the deviation angle is positive when the blade should push, and negative when the blade should pull
    FloatType deviationAngle = -std::asin(normalVelocity / relativeVelocity.magnitude());
    return deviationAngle;
}

Vector3D CSourceBFMBase::computeInviscidForceDirection(const Vector3D& relativeVelocity, const Vector3D& normalCamber){
    Vector3D relativeVelocityPlus = {relativeVelocity.x(), relativeVelocity.y(), relativeVelocity.z()};
    Vector3D wDir = relativeVelocityPlus.normalized();
    Vector3D normal = normalCamber.normalized();
    
    // method based on versor construction, used as a first version
    // Vector3D relativeVelocityPlus = {relativeVelocity.x() + 1E-6, relativeVelocity.y() + 1E-6, relativeVelocity.z() + 1E-6};
    // FloatType A, B, C, Delta;
    // A = wDir.z()*wDir.z() + wDir.x()*wDir.x();
    // B = 2 * wDir.y() * wDir.x() * normal.y();
    // C = (wDir.z()*wDir.z() * normal.y()*normal.y()) + (wDir.y()*wDir.y() * normal.y()*normal.y()) - wDir.z()*wDir.z();
    // Delta = B*B - 4*A*C;

    // if (Delta < 0) {
    //     std::cout << "Inviscid force direction computation found no real result. Radial component set to zero." << std::endl;
    //     Vector3D versor = {-wDir.z(), 0, wDir.x()};
    //     return versor.normalized();
    // }

    // FloatType fAxial1 = (-B + std::sqrt(Delta)) / (2 * A);
    // FloatType fAxial2 = (-B - std::sqrt(Delta)) / (2 * A);

    // FloatType fTangential1 = computeTangentialComponent(fAxial1, wDir, normal);
    // FloatType fTangential2 = computeTangentialComponent(fAxial2, wDir, normal);

    // FloatType fRadial1 = normal.y();
    // FloatType fRadial2 = normal.y();

    // Vector3D versor1 = {fAxial1, fRadial1, fTangential1};
    // Vector3D versor2 = {fAxial2, fRadial2, fTangential2};

    // if (versor1.dot(normal) > 0) {
    //     return versor1;
    // } else {
    //     return versor2;
    // }

    // method based on versor projection
    // the idea is that the loading versor is found be subtracting from the normal camber versor its projection on the relative velocity versor
    Vector3D versor = normal - wDir * normal.dot(wDir);
    versor = versor.normalized();

    return versor;

}


Vector3D CSourceBFMBase::computeInviscidForceDirection(const Vector3D& relativeVelocity, const FloatType& gaspathAngle, const FloatType& flowAngle, FloatType& leanAngle){
    // overwrite lean to zero for now
    // leanAngle = 0.0;
    
    // the minus are there because the derivation considers opposite signs for the rotations of the ref frames
    FloatType nAxial = -std::sin(flowAngle) * std::cos(leanAngle) * std::cos(gaspathAngle) - std::sin(leanAngle) * std::sin(gaspathAngle);
    FloatType nRadial = -std::sin(flowAngle) * std::sin(gaspathAngle) * std::cos(leanAngle) + std::sin(leanAngle) * std::cos(gaspathAngle);
    FloatType nTan = std::cos(flowAngle) * std::cos(leanAngle);

    Vector3D versor =   {nAxial, nRadial, nTan};

    // convention on the sign for a positive pushing force
    if (_flowAngle > 0){
        versor = -versor;
    }

    return versor;
}


FloatType CSourceBFMBase::computeTangentialComponent(FloatType fAxial, const Vector3D& relativeVelocityDirection, const Vector3D& normalCamber){
    FloatType fTangential = (-relativeVelocityDirection.x() * fAxial - relativeVelocityDirection.y() * normalCamber.y()) / relativeVelocityDirection.z();
    return fTangential;
}


FloatType CSourceBFMBase::computeBodyForcePerturbationScaling(size_t i, size_t j, size_t k, FloatType timePhysical) const{

    if (!_config.isPerturbationBodyForceActive()){
        return 0.0;
    }

    // if time not in the perturbation, return 1
    if (timePhysical < _perturbationTimeStart || timePhysical > _perturbationTimeStart + _perturbationTimeDuration){
        return 0.0;
    }

    // if the rest is not true, compute the scaling factor of the body force, through RBF in space and time
    FloatType variationMax = _perturbationScalingFactor - 1.0;

    Vector3D point = _mesh.getVertex(i,j,k);
    FloatType d = (point - _perturbationCenter).magnitude(); // distance of the point from the perturbation center

    if (d > 2.0*_perturbationRadius){
        return 0.0;
    }

    FloatType spaceModulation = std::exp(-(d*d)/(_perturbationRadius*_perturbationRadius)); // RBF in space. At d=R the modulation is around 35%

    FloatType timeCentral = _perturbationTimeStart + _perturbationTimeDuration/3.0; // central time touches the peak
    FloatType sigma = _perturbationTimeDuration/2.0; // sharp rise and fall
    FloatType timeModulation = std::exp( - (std::pow(timePhysical-timeCentral, 2)) / std::pow(sigma,2)); // RBF in time

    FloatType perturbationVariation = variationMax * spaceModulation * timeModulation;

    // if (spaceModulation > 0.5){
    //     std::cout<< "For point (" << i << "," << j << "," << k << "):" << std::endl;
    //     std::cout << "Perturbation space scaling factor: " << spaceModulation << std::endl;
    // }

    // if (timeModulation > 0.5){
    //     std::cout<< "For point (" << i << "," << j << "," << k << "):" << std::endl;
    //     std::cout << "Perturbation time scaling factor: " << timeModulation << std::endl;
    // }

    return perturbationVariation;
}
