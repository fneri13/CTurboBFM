#include "CSourceBFMBase.hpp"


CSourceBFMBase::CSourceBFMBase(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : _config(config), _fluid(fluid), _mesh(mesh) {

            // read data associated to body force perturbation
            _perturbationBodyForceActive = _config.isPerturbationBodyForceActive();
            if (_perturbationBodyForceActive){
                
                _perturbationCenterI = _config.getPerturbationIJK_Coords()[0];
                _perturbationCenterJ = _config.getPerturbationIJK_Coords()[1];
                _perturbationCenterK = _config.getPerturbationIJK_Coords()[2];

                if (_perturbationCenterI > _mesh.getNumberPointsI()-1 || _perturbationCenterJ > _mesh.getNumberPointsJ()-1 || _perturbationCenterK > _mesh.getNumberPointsK()-1){
                    throw std::runtime_error("Perturbation center is outside the mesh.");
                }

                _perturbationExtensionI = _config.getPerturbationIJK_Extension()[0];
                _perturbationExtensionJ = _config.getPerturbationIJK_Extension()[1];
                _perturbationExtensionK = _config.getPerturbationIJK_Extension()[2];

                if (_perturbationCenterI+_perturbationExtensionI > _mesh.getNumberPointsI()-1 ||
                    _perturbationCenterJ+_perturbationExtensionJ > _mesh.getNumberPointsJ()-1 || 
                    _perturbationCenterK+_perturbationExtensionK > _mesh.getNumberPointsK()-1){
                    
                        throw std::runtime_error("Perturbation extension is outside the mesh.");
                }

                if (_perturbationCenterI-_perturbationExtensionI > _mesh.getNumberPointsI()-1 || 
                    _perturbationCenterJ-_perturbationExtensionJ > _mesh.getNumberPointsJ()-1 || 
                    _perturbationCenterK-_perturbationExtensionK > _mesh.getNumberPointsK()-1){
                    
                        throw std::runtime_error("Perturbation extension is outside the mesh.");
                }
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
    FloatType bodyForceScaling = computeBodyForcePerturbationScaling(i, j, k, timePhysical);

    deviationAngle(i, j, k) = _deviationAngle; // update the deviation angle that has been computed from every bfm model, getting it ready for the output

    return (blockageSource + bodyForceSource) * bodyForceScaling;
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
    
    // blade properties
    _numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    _pitch = 2.0 * M_PI * _radius / _numberBlades;
    _normalCamberAxial = _mesh.getInputFields(FieldNames::NORMAL_AXIAL, i, j, k);
    _normalCamberRadial = _mesh.getInputFields(FieldNames::NORMAL_RADIAL, i, j, k);
    _normalCamberTangential = _mesh.getInputFields(FieldNames::NORMAL_TANGENTIAL, i, j, k);
    _normalCamberCylindric = {_normalCamberAxial, _normalCamberRadial, _normalCamberTangential};
    _blockage = _mesh.getInputFields(FieldNames::BLOCKAGE, i, j, k);
    _bladeIsPresent = _mesh.getInputFields(FieldNames::BLADE_PRESENT, i, j, k);
    
    // force directions 
    _deviationAngle = computeDeviationAngle(_relativeVelocityCylindric, _normalCamberCylindric);
    _inviscidForceDirectionCylindrical = computeInviscidForceDirection(_relativeVelocityCylindric, _normalCamberCylindric);
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
    Vector3D relativeVelocityPlus = {relativeVelocity.x() + 1E-6, relativeVelocity.y() + 1E-6, relativeVelocity.z() + 1E-6};
    Vector3D wDir = relativeVelocityPlus.normalized();
    Vector3D normal = normalCamber.normalized();
    
    FloatType A, B, C, Delta;
    A = wDir.z()*wDir.z() + wDir.x()*wDir.x();
    B = 2 * wDir.y() * wDir.x() * normal.y();
    C = (wDir.z()*wDir.z() * normal.y()*normal.y()) + (wDir.y()*wDir.y() * normal.y()*normal.y()) - wDir.z()*wDir.z();
    Delta = B*B - 4*A*C;

    if (Delta < 0) {
        std::cout << "Inviscid force direction computation found no real result. Radial component set to zero." << std::endl;
        Vector3D versor = {-wDir.z(), 0, wDir.x()};
        return versor.normalized();
    }

    FloatType fAxial1 = (-B + std::sqrt(Delta)) / (2 * A);
    FloatType fAxial2 = (-B - std::sqrt(Delta)) / (2 * A);

    FloatType fTangential1 = computeTangentialComponent(fAxial1, wDir, normal);
    FloatType fTangential2 = computeTangentialComponent(fAxial2, wDir, normal);

    FloatType fRadial1 = normal.y();
    FloatType fRadial2 = normal.y();

    Vector3D versor1 = {fAxial1, fRadial1, fTangential1};
    Vector3D versor2 = {fAxial2, fRadial2, fTangential2};

    if (versor1.dot(normal) > 0) {
        return versor1;
    } else {
        return versor2;
    }
}


FloatType CSourceBFMBase::computeTangentialComponent(FloatType fAxial, const Vector3D& relativeVelocityDirection, const Vector3D& normalCamber){
    FloatType fTangential = (-relativeVelocityDirection.x() * fAxial - relativeVelocityDirection.y() * normalCamber.y()) / relativeVelocityDirection.z();
    return fTangential;
}


FloatType CSourceBFMBase::computeBodyForcePerturbationScaling(size_t i, size_t j, size_t k, FloatType timePhysical) const{

    // if time not in the perturbation, return 1
    if (timePhysical < _perturbationTimeStart || timePhysical > _perturbationTimeStart + _perturbationTimeDuration){
        return 1.0;
    }

    // if location out of perturbation, return 1
    if (i > _perturbationCenterI + _perturbationExtensionI || 
        j > _perturbationCenterJ + _perturbationExtensionJ || 
        k > _perturbationCenterK + _perturbationExtensionK ||
        i < _perturbationCenterI - _perturbationExtensionI || 
        j < _perturbationCenterJ - _perturbationExtensionJ || 
        k < _perturbationCenterK - _perturbationExtensionK){
        
            return 1.0;
    }

    return _perturbationScalingFactor;
}
