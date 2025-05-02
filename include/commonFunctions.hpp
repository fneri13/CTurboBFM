#pragma once

#include "types.hpp"
#include "CFluid.hpp"

// compute the integral flux of an array2D of vectors through an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

// compute the total area of an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces);

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<FloatType> &vecX, const Matrix2D<FloatType> &vecY, const Matrix2D<FloatType> &vecZ);


StateVector getEulerPrimitiveFromConservative(StateVector conservative);

StateVector getEulerConservativeFromPrimitive(StateVector primitive);

StateVector computeEulerFluxFromPrimitive(StateVector primitive, Vector3D surface, CFluid fluid);

StateVector computeEulerFluxFromConservative(StateVector conservative, Vector3D surface, CFluid fluid);

FloatType computeAngleBetweenVectors(const Vector3D& v1, const Vector3D& v2);

Vector3D rotateVectorAlongXAxis(const Vector3D& v, FloatType theta);

// convert a cartesian vector (x,y,z components) to a cylindrical vector (axial, radial, tangential components)
Vector3D computeCylindricalVectorFromCartesian(Vector3D vec, FloatType theta);

// convert a cylindrical vector (axial, radial, tangential components) to a cartesian vector (x,y,z components)
Vector3D computeCartesianVectorFromCylindrical(Vector3D vec, FloatType theta);

// integrate the radial equilibrium equation using trapezoidal rule
void integrateRadialEquilibrium(const std::vector<FloatType>& density, const std::vector<FloatType>& velocityTang, 
                                const std::vector<FloatType>& radius, const FloatType& hubPressure, 
                                std::vector<FloatType>& pressure);
