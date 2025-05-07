#include "CSourceBFMBase.hpp"

StateVector CSourceBFMBase::computeSource(size_t i, size_t j, size_t k, const StateVector& primitive) const {
    StateVector blockageSource = computeBlockageSource(i, j, k, primitive);
    StateVector bodyForceSource = computeBodyForceSource(i, j, k, primitive);
    return blockageSource + bodyForceSource;
}

StateVector CSourceBFMBase::computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive) const {
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


StateVector CSourceBFMBase::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) const {
    return StateVector({0, 0, 0, 0, 0});
}