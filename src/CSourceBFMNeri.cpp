#include "CSourceBFMNeri.hpp"

StateVector CSourceBFMNeri::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMNeri::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType fn_0 = _mesh.getInputFields(FieldNames::INFERENCE_FN_0, i, j, k);
    FloatType fn_1 = _mesh.getInputFields(FieldNames::INFERENCE_FN_1, i, j, k);
    FloatType fn_2 = _mesh.getInputFields(FieldNames::INFERENCE_FN_2, i, j, k);
    FloatType fn_3 = _mesh.getInputFields(FieldNames::INFERENCE_FN_3, i, j, k);
    
    FloatType pressure = _fluid.computePressure_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);

    // polynomial inference
    FloatType forceMag = fn_0 + (fn_1 * pressure) + (fn_2 * pressure * pressure) + (fn_3 * pressure * pressure * pressure);

    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}



StateVector CSourceBFMNeri::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType fp_0 = _mesh.getInputFields(FieldNames::INFERENCE_FP_0, i, j, k);
    FloatType fp_1 = _mesh.getInputFields(FieldNames::INFERENCE_FP_1, i, j, k);
    FloatType fp_2 = _mesh.getInputFields(FieldNames::INFERENCE_FP_2, i, j, k);
    FloatType fp_3 = _mesh.getInputFields(FieldNames::INFERENCE_FP_3, i, j, k);
    
    FloatType pressure = _fluid.computePressure_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);

    // polynomial inference
    FloatType forceMag = fp_0 + (fp_1 * pressure) + (fp_2 * pressure * pressure) + (fp_3 * pressure * pressure * pressure);

    Vector3D forceCylindrical = _viscousForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}


