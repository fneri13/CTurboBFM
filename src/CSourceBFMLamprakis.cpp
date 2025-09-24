#include "CSourceBFMLamprakis.hpp"

CSourceBFMLamprakis::CSourceBFMLamprakis(const Config &config, const CFluid &fluid, const CMesh &mesh,
                                  std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance) 
    : CSourceBFMBase(config, fluid, mesh), _turboPerformance(turboPerformance)
{}


StateVector CSourceBFMLamprakis::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);

    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}

void CSourceBFMLamprakis::computeFlowState(size_t i, size_t j, size_t k, const StateVector& primitive) {
    // geometrical info
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
    _relativeVelocityMag = _relativeVelocityCartesian.magnitude();
    
    // blade properties
    _numberBlades = _mesh.getInputFields(FieldNames::NUMBER_BLADES, i, j, k);
    _blockage = _mesh.getInputFields(FieldNames::BLOCKAGE, i, j, k);
    _bladeIsPresent = _mesh.getInputFields(FieldNames::BLADE_PRESENT, i, j, k);
    _metalAngle = _mesh.getInputFields(FieldNames::BLADE_METAL_ANGLE, i, j, k);
    _d_metalAngle_dm = _mesh.getInputFields(FieldNames::D_BLADE_METAL_ANGLE_DM, i, j, k);
    _leanAngle = _mesh.getInputFields(FieldNames::BLADE_LEAN_ANGLE, i, j, k);
    _gasPathAngle = _mesh.getInputFields(FieldNames::BLADE_GAS_PATH_ANGLE, i, j, k);

    // flow kinematics and angles
    _velMeridionalMag = std::sqrt(_relativeVelocityCylindric.x() * _relativeVelocityCylindric.x() + _relativeVelocityCylindric.y() * _relativeVelocityCylindric.y());
    _betaFlow = std::atan2(_relativeVelocityCylindric.z(), _velMeridionalMag);

    // for the moment the flow is perfectly aligned with the camber line
    _loadingVersor_cyl = computeLoadingVersor(_metalAngle, _leanAngle, _gasPathAngle);
    _loadingVersor_cart = computeCartesianVectorFromCylindrical(_loadingVersor_cyl, _theta);

}

Vector3D CSourceBFMLamprakis::computeLoadingVersor(const FloatType& metalAngle, const FloatType& leanAngle, const FloatType& gasPathAngle) const {
    // Later the metal angle needs to be corrected with the deviation (or slip)
    Vector3D versor {0,0,0};
    versor.x() = -std::sin(-metalAngle) * std::cos(gasPathAngle);
    versor.y() = std::sin(-metalAngle) * std::sin(gasPathAngle);
    versor.z() = std::cos(-metalAngle);
    return versor;
}


StateVector CSourceBFMLamprakis::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType fn = computeLoadingForceIntensity(primitive);
    inviscidForce(i, j, k) = _loadingVersor_cart * fn * _bladeIsPresent;
    Vector3D forceCylindrical = computeCylindricalVectorFromCartesian(inviscidForce(i, j, k), _theta);
    
    StateVector source({0,0,0,0,0});
    source[1] = inviscidForce(i, j, k).x();
    source[2] = inviscidForce(i, j, k).y();
    source[3] = inviscidForce(i, j, k).z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

FloatType CSourceBFMLamprakis::computeLoadingForceIntensity(const StateVector& primitive) const {
    FloatType fn = 0;
    FloatType w = _relativeVelocityMag;
    FloatType soundSpeed = _fluid.computeSoundSpeed_rho_u_et(primitive[0], _relativeVelocityCartesian, primitive[4]);
    FloatType Mrel = _velMeridionalMag / soundSpeed;
    FloatType Mom = _omega*_radius / soundSpeed;
    FloatType gmma = _fluid.getGamma();
    FloatType beta = _metalAngle;
    FloatType dbdm = _d_metalAngle_dm;
    FloatType phi = _gasPathAngle;

    FloatType term1, term2, term3, term4;
    beta *= -1;
    term1 = w*w * (std::cos(beta) * ((2.0*gmma-1.0)*Mrel*Mrel - std::tan(beta)*std::tan(beta) + 1.0) * dbdm / ((2*gmma-1.0)*Mrel*Mrel + 1.0) + 
                (Mom*Mom/Mrel/Mrel + (2.0*gmma-1.0)*Mrel*Mrel + 2.0) * std::sin(phi)*std::sin(beta) / _radius / ((2*gmma-1.0)*Mrel*Mrel + 1.0));
    
    term2 = 2.0 * w * _omega * std::sin(phi);

    term3 = _omega*_omega * _radius * std::sin(phi) * std::sin(beta);

    term4 = 0.0; //for now

    fn = term1 + term2 + term3 + term4;
    return fn;
}


StateVector CSourceBFMLamprakis::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    viscousForce(i, j, k) = {0, 0, 0};
    Vector3D forceCylindrical = computeCylindricalVectorFromCartesian(viscousForce(i, j, k), _theta);
    
    StateVector source({0,0,0,0,0});
    source[1] = viscousForce(i, j, k).x();
    source[2] = viscousForce(i, j, k).y();
    source[3] = viscousForce(i, j, k).z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}
    



