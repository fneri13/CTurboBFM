#include "CSourceBFMBase.hpp"

StateVector CSourceBFMBase::computeSource(size_t i, size_t j, size_t k){
    StateVector blockageSource = computeBlockageSource(i, j, k);
    StateVector bodyForceSource = computeBodyForceSource(i, j, k);
    return blockageSource + bodyForceSource;
}


StateVector CSourceBFMBase::computeBlockageSource(size_t i, size_t j, size_t k){
    StateVector conservative = _solution.at(i,j,k);
    StateVector primitive = getEulerPrimitiveFromConservative(conservative);

    Vector3D velocity({primitive[1], primitive[2], primitive[3]});
    FloatType totalEnthalpy = _fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType blockage = _mesh.getInputFields(FieldNames::BLOCKAGE, i, j, k);
    return StateVector({0, 0, 0, -blockage * totalEnthalpy});

}

StateVector CSourceBFMBase::computeBodyForceSource(size_t i, size_t j, size_t k){
    return StateVector({0, 0, 0, 0, 0});
}