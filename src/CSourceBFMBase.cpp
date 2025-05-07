#include "CSourceBFMBase.hpp"

StateVector CSourceBFMBase::computeSource(size_t i, size_t j, size_t k, const StateVector& primitive) {

    StateVector blockageSource({0,0,0,0,0});
    if (_config.isBlockageActive()){
        blockageSource = computeBlockageSource(i, j, k, primitive);
    }

    FloatType numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    if (numberBlades == 0){ // this is the layer upstream of leading edge, no blade present, only blockage contribution
        return blockageSource;
    }

    StateVector bodyForceSource = computeBodyForceSource(i, j, k, primitive);
    return blockageSource + bodyForceSource;
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


StateVector CSourceBFMBase::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) {
    return StateVector({0, 0, 0, 0, 0});
}

void CSourceBFMBase::computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive){
    _point = _mesh.getVertex(i, j, k);
    _radius = std::sqrt(_point.z() * _point.z() + _point.y() * _point.y());
    _theta = std::atan2(_point.z(), _point.y());
    _velocityCartesian = {primitive[1], primitive[2], primitive[3]};
    _velocityCylindrical = computeCylindricalVectorFromCartesian(_velocityCartesian, _theta);
    _omega = _mesh.getInputFields(FieldNames::RPM, i, j, k) * 2 * M_PI / 60;
    _dragVelocityCylindrical = {0, 0, _omega * _radius};
    _relativeVelocityCylindric = _velocityCylindrical - _dragVelocityCylindrical;
    _numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    _pitch = 2.0 * M_PI * _radius / _numberBlades;
    _normalCamberAxial = _mesh.getInputFields(FieldNames::NORMAL_AXIAL, i, j, k);
    _normalCamberRadial = _mesh.getInputFields(FieldNames::NORMAL_RADIAL, i, j, k);
    _normalCamberTangential = _mesh.getInputFields(FieldNames::NORMAL_TANGENTIAL, i, j, k);
    _normalCamberCylindric = {_normalCamberAxial, _normalCamberRadial, _normalCamberTangential};
    _deviationAngle = computeDeviationAngle(_relativeVelocityCylindric, _normalCamberCylindric);
    _inviscidForceDirection = computeInviscidForceDirection(_relativeVelocityCylindric, _normalCamberCylindric);
    _relativeVelocityCartesian = computeCartesianVectorFromCylindrical(_relativeVelocityCylindric, _theta);
    _viscousForceDirection = - _relativeVelocityCartesian / _relativeVelocityCartesian.magnitude();
}


FloatType CSourceBFMBase::computeDeviationAngle(Vector3D relativeVelocity, Vector3D normalCamber){
    Vector3D normal = normalCamber / normalCamber.magnitude();
    FloatType normalVelocity = relativeVelocity.dot(normal);
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
