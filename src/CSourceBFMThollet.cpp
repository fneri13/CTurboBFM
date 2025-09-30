#include "CSourceBFMThollet.hpp"

StateVector CSourceBFMThollet::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMThollet::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType Kmach = computeCompressibilityCorrection(_relativeVelocityCylindric, primitive);
    FloatType forceMag = _Kn * Kmach * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * M_PI * _deviationAngle / _pitch / std::abs(_normalCamberTangential) / _blockage;
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

FloatType CSourceBFMThollet::computeCompressibilityCorrection(const Vector3D& relativeVelocityCylindric, const StateVector& primitive) {
    FloatType soundSpeed = _fluid.computeSoundSpeed_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
    FloatType relativeMach = relativeVelocityCylindric.magnitude() / soundSpeed;

    if (std::abs(relativeMach-1.0) < 0.001)
        relativeMach = 0.999;
    
    FloatType kPrime = 0.0;
    
    if (relativeMach<1.0)
        kPrime = 1.0 / std::sqrt(1.0 - relativeMach*relativeMach);
    else
        kPrime = 2.0 / M_PI / std::sqrt(relativeMach*relativeMach - 1.0);

    return std::min(3.0, kPrime);
}

StateVector CSourceBFMThollet::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType nu = _config.getFluidKinematicViscosity();
    FloatType stwl = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i, j, k);

    // if the stwl is zero, get the average with the downstream points, otherwise the ReX will be zero, and therefore Cf infinity
    if (stwl < 1e-9){
        stwl = (_mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i+1, j, k) + _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i, j, k)) /2.0 ;
    }

    FloatType ReX = _relativeVelocityCylindric.magnitude() * stwl / nu;
    FloatType Cf = 0.0592 * std::pow(ReX, -0.2);

    FloatType forceMag = 0.0;
    if (_offDesignActive == "no"){
        forceMag = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / (_pitch * _blockage * std::abs(_normalCamberTangential)) * (_Kf * Cf);
        }
    else if (_offDesignActive == "yes"){
        FloatType deviationAnglePivot = _mesh.getInputFields(FieldNames::DEVIATION_ANGLE_PIVOT, i, j, k);
        FloatType Kmach = computeCompressibilityCorrection(_relativeVelocityCylindric, primitive);

        // FloatType deltaDev = std::abs(_deviationAngle - deviationAnglePivot);
        // FloatType term1 = _Kf * Cf;
        // FloatType term2 = M_PI * Kmach * std::pow(std::abs(_deviationAngle - deviationAnglePivot), 2.0 * _Kd);

        forceMag = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / (_pitch * _blockage * std::abs(_normalCamberTangential)) * (
            _Kf * Cf + M_PI * Kmach * std::pow(std::abs(_deviationAngle - deviationAnglePivot), 2.0 * _Kd));
    }
    else {
        throw std::runtime_error("HALL_THOLLET_OFF_DESIGN_ACTIVE can only be <yes> or <no>");
    }

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