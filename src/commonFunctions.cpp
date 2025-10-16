#include "commonFunctions.hpp"
#include <iostream>


FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).dot(vectors(i,j));
        }
    }
    return sum;
}

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<FloatType> &vecX, const Matrix2D<FloatType> &vecY, const Matrix2D<FloatType> &vecZ){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).x() * vecX(i,j) + surfaces(i,j).y() * vecY(i,j) + surfaces(i,j).z() * vecZ(i,j);
        }
    }
    return sum;
}


FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (size_t i=0; i<ni; i++){
        for (size_t j=0; j<nj; j++){
            sum += surfaces(i,j).magnitude();
        }
    }
    return sum;
}


StateVector getEulerPrimitiveFromConservative(StateVector conservative){
    StateVector primitive;
    primitive[0] = conservative[0];
    primitive[1] = conservative[1] / conservative[0];
    primitive[2] = conservative[2] / conservative[0];
    primitive[3] = conservative[3] / conservative[0];
    primitive[4] = conservative[4] / conservative[0];
    return primitive;
}


StateVector getEulerConservativeFromPrimitive(StateVector primitive){
    StateVector conservative;
    conservative[0] = primitive[0];
    conservative[1] = primitive[1] * primitive[0];
    conservative[2] = primitive[2] * primitive[0];
    conservative[3] = primitive[3] * primitive[0];
    conservative[4] = primitive[4] * primitive[0];
    return conservative;
}


StateVector computeEulerFluxFromPrimitive(StateVector primitive, Vector3D surface, const CFluidBase& fluid){
    Vector3D normal = surface / surface.magnitude();
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    FloatType normalVelocity = velocity.dot(normal);
    FloatType pressure = fluid.computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType totEnthalpy = fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);

    StateVector flux {};
    flux[0] = primitive[0] * normalVelocity;
    flux[1] = primitive[0] * normalVelocity * velocity.x() + pressure * normal.x();
    flux[2] = primitive[0] * normalVelocity * velocity.y() + pressure * normal.y();
    flux[3] = primitive[0] * normalVelocity * velocity.z() + pressure * normal.z();
    flux[4] = primitive[0] * normalVelocity * totEnthalpy;
    return flux;

}


StateVector computeEulerFluxFromConservative(StateVector conservative, Vector3D surface, const CFluidBase& fluid){
    StateVector primitive = getEulerPrimitiveFromConservative(conservative);
    StateVector flux = computeEulerFluxFromPrimitive(primitive, surface, fluid);
    return flux;
}


FloatType computeAngleBetweenVectors(const Vector3D& v1, const Vector3D& v2) {
    Vector3D v1norm = v1.normalized();
    Vector3D v2norm = v2.normalized();

    FloatType dot = v1norm.dot(v2norm);

    // Clamp dot product to valid range for acos()
    dot = std::clamp(dot, static_cast<FloatType>(-1.0), static_cast<FloatType>(1.0));

    return std::acos(dot);
}


Vector3D rotateVectorAlongXAxis(const Vector3D& v, FloatType theta){
    Vector3D newVector;
    newVector.x() = v.x();
    newVector.y() = v.y() * cos(theta) - v.z() * sin(theta);
    newVector.z() = v.y() * sin(theta) + v.z() * cos(theta);
    return newVector;
}

Vector3D computeCylindricalVectorFromCartesian(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y() * cos(theta) + vec.z() * sin(theta);
    newVector.z() = - vec.y() * sin(theta) + vec.z() * cos(theta);
    return newVector;
}

Vector3D computeCartesianVectorFromCylindrical(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y() * cos(theta) - vec.z() * sin(theta);
    newVector.z() = vec.y() * sin(theta) + vec.z() * cos(theta);
    return newVector;
}


void integrateRadialEquilibrium(const std::vector<FloatType>& density,
    const std::vector<FloatType>& velocityTang,
    const std::vector<FloatType>& radius,
    const FloatType& hubPressure,
    std::vector<FloatType>& pressure)
{

    const size_t N = radius.size();
    pressure[0] = hubPressure;

    for (size_t i = 1; i < N; ++i) {
        FloatType r1 = radius[i - 1];
        FloatType r2 = radius[i];

        FloatType dr = r2 - r1;

        FloatType rho1 = density[i - 1];
        FloatType rho2 = density[i];

        FloatType ut1 = velocityTang[i - 1];
        FloatType ut2 = velocityTang[i];

        // Evaluate dp/dr = rho * u_theta^2 / r at both points
        FloatType dpdr1 = rho1 * ut1 * ut1 / r1;
        FloatType dpdr2 = rho2 * ut2 * ut2 / r2;

        // Trapezoidal integration step
        pressure[i] = pressure[i - 1] + 0.5 * (dpdr1 + dpdr2) * dr;
    }

}


void computeGradientGreenGauss(const Matrix3D<Vector3D>& surfacesI, const Matrix3D<Vector3D>& surfacesJ, const Matrix3D<Vector3D>& surfacesK, 
    const Matrix3D<Vector3D>& midPointsI, const Matrix3D<Vector3D>& midPointsJ, const Matrix3D<Vector3D>& midPointsK, 
    const Matrix3D<Vector3D>& nodes, const Matrix3D<
    FloatType>& volumes, Matrix3D<FloatType>& field, Matrix3D<Vector3D>& gradient){

    size_t ni = volumes.sizeI();
    size_t nj = volumes.sizeJ();
    size_t nk = volumes.sizeK();

    for (size_t i = 0; i < ni; i++){
        for (size_t j = 0; j < nj; j++){
            for (size_t k = 0; k < nk; k++){
                Vector3D surfaceWest = - surfacesI(i,j,k);
                Vector3D surfaceEast = surfacesI(i+1,j,k);
                Vector3D surfaceNorth = surfacesJ(i,j+1,k);
                Vector3D surfaceSouth = - surfacesJ(i,j,k);
                Vector3D surfaceBottom = - surfacesK(i,j,k);
                Vector3D surfaceTop = surfacesK(i,j,k+1);

                Vector3D centerWest = midPointsI(i,j,k);
                Vector3D centerEast = midPointsI(i+1,j,k);
                Vector3D centerNorth = midPointsJ(i,j+1,k);
                Vector3D centerSouth = midPointsJ(i,j,k);
                Vector3D centerBottom = midPointsK(i,j,k);
                Vector3D centerTop = midPointsK(i,j,k+1);

                FloatType scalarWest = interpolateScalar(field, i, j, k, nodes, centerWest, Direction3D::WEST);
                FloatType scalarEast = interpolateScalar(field, i, j, k, nodes, centerEast, Direction3D::EAST);
                FloatType scalarNorth = interpolateScalar(field, i, j, k, nodes, centerNorth, Direction3D::NORTH);
                FloatType scalarSouth = interpolateScalar(field, i, j, k, nodes, centerSouth, Direction3D::SOUTH);
                FloatType scalarBottom = interpolateScalar(field, i, j, k, nodes, centerBottom, Direction3D::BOTTOM);
                FloatType scalarTop = interpolateScalar(field, i, j, k, nodes, centerTop, Direction3D::TOP);

                std::array<Vector3D, 6> surfaces = {surfaceWest, surfaceEast, surfaceNorth, surfaceSouth, surfaceBottom, surfaceTop};
                std::array<FloatType, 6> scalars = {scalarWest, scalarEast, scalarNorth, scalarSouth, scalarBottom, scalarTop};

                gradient(i,j,k) = computeGreenGaussFormula(surfaces, scalars, volumes(i,j,k));
            }
        }
    }
}


FloatType interpolateScalar(const Matrix3D<FloatType>& field, size_t i, size_t j, size_t k, const Matrix3D<Vector3D>& nodes, const Vector3D& point, Direction3D direction){
    Vector3D point0 = nodes(i,j,k);
    FloatType field0 = field(i,j,k);
    
    FloatType field1 = 0.0;
    Vector3D point1(0.0,0.0,0.0);
    
    FloatType fieldInterpolated = 0.0;

    size_t ni = field.sizeI();
    size_t nj = field.sizeJ();
    size_t nk = field.sizeK();

    if (direction==Direction3D::WEST && i==0){
        return field0;
    }
    else if (direction==Direction3D::WEST){
        point1 = nodes(i-1,j,k);
        field1 = field(i-1,j,k);
    }
    else if (direction==Direction3D::EAST && i==ni-1){
        return field0;
    }
    else if (direction==Direction3D::EAST){
        point1 = nodes(i+1,j,k);
        field1 = field(i+1,j,k);
    }

    else if (direction==Direction3D::NORTH && j==nj-1){
        return field0;
    }
    else if (direction==Direction3D::NORTH){
        point1 = nodes(i,j+1,k);
        field1 = field(i,j+1,k);
    }
    else if (direction==Direction3D::SOUTH && j==0){
        return field0;
    }
    else if (direction==Direction3D::SOUTH){
        point1 = nodes(i,j-1,k);
        field1 = field(i,j-1,k);
    }

    else if (direction==Direction3D::BOTTOM && k==0){
        return field0;
    }
    else if (direction==Direction3D::BOTTOM){
        point1 = nodes(i,j,k-1);
        field1 = field(i,j,k-1);
    }
    else if (direction==Direction3D::TOP && k==nk-1){
        return field0;
    }
    else if (direction==Direction3D::TOP){
        point1 = nodes(i,j,k+1);
        field1 = field(i,j,k+1);
    }

    FloatType distance1 = (point - point0).magnitude();
    FloatType distance2 = (point - point1).magnitude();
    fieldInterpolated = field0 + (field1 - field0) * distance1 / (distance1 + distance2);
    return fieldInterpolated;
}


Vector3D computeGreenGaussFormula(const std::array<Vector3D,6>& surfaces, const std::array<FloatType,6>& scalars, const FloatType& volume){
    Vector3D gradient(0.0,0.0,0.0);
    int length = surfaces.size();
    for (int i=0; i<length; i++){
        gradient += surfaces[i] * scalars[i];
    }
    return gradient / volume;
}

FloatType atan2_from0_to2pi(FloatType y, FloatType x){
    FloatType angle = std::atan2(y, x);
    if (angle < 0)
        angle += 2 * M_PI;
    return angle;
}


FloatType inverseRescalingMinMax(const FloatType& value, const FloatType& min, const FloatType& max){
    FloatType res = min + (max - min) * value;
    return res;
}


FloatType rescaleMinMax(const FloatType& value, const FloatType& min, const FloatType& max){
    FloatType res = (value - min) / (max - min);
    return res;
}




FloatType interpolateLinear(const std::vector<double>& x,
                            const std::vector<double>& y,
                            const FloatType& xp) {

    // Extrapolate to the left
    if (xp <= x.front()) {
        double slope = (y[1] - y[0]) / (x[1] - x[0]);
        return y[0] + slope * (xp - x[0]);
    }

    // Extrapolate to the right
    if (xp >= x.back()) {
        size_t n = x.size();
        double slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
        return y[n - 1] + slope * (xp - x[n - 1]);
    }

    // Find the interval [x[i], x[i+1]] such that x[i] <= xp < x[i+1]
    auto it = std::upper_bound(x.begin(), x.end(), xp);
    size_t i = std::distance(x.begin(), it) - 1;

    double x0 = x[i];
    double x1 = x[i + 1];
    double y0 = y[i];
    double y1 = y[i + 1];

    double t = (xp - x0) / (x1 - x0);
    return ((1.0 - t) * y0 + t * y1);
}
