#pragma once

#include "types.hpp"
#include "CFluidBase.hpp"
#include "CFluidIdeal.hpp"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>

// compute the integral flux of an array2D of vectors through an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

// compute the total area of an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces);

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<FloatType> &vecX, const Matrix2D<FloatType> &vecY, const Matrix2D<FloatType> &vecZ);


StateVector getEulerPrimitiveFromConservative(StateVector conservative);

StateVector getEulerConservativeFromPrimitive(StateVector primitive);

StateVector computeEulerFluxFromPrimitive(StateVector primitive, Vector3D surface, const CFluidBase& fluid);

StateVector computeEulerFluxFromConservative(StateVector conservative, Vector3D surface, const CFluidBase& fluid);

Matrix2D<FloatType> computeAdvectionFluxJacobian(const StateVector& primitive, const Vector3D& direction, const CFluidBase& fluid);

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


void computeGradientGreenGauss(const Matrix3D<Vector3D>& surfacesI, const Matrix3D<Vector3D>& surfacesJ, const Matrix3D<Vector3D>& surfacesK, 
                                const Matrix3D<Vector3D>& midPointsI, const Matrix3D<Vector3D>& midPointsJ, const Matrix3D<Vector3D>& midPointsK, 
                                const Matrix3D<Vector3D>& nodes, const Matrix3D<FloatType>& volumes, const Matrix3D<FloatType>& field, Matrix3D<Vector3D>& gradient);


FloatType interpolateScalar(const Matrix3D<FloatType>& field, size_t i, size_t j, size_t k, const Matrix3D<Vector3D>& nodes, const Vector3D& point, Direction3D direction);

Vector3D computeGreenGaussFormula(const std::array<Vector3D,6>& surfaces, const std::array<FloatType,6>& scalars, const FloatType& volume);

FloatType atan2_from0_to2pi(FloatType y, FloatType x);

FloatType inverseRescalingMinMax(const FloatType& value, const FloatType& min, const FloatType& max);

FloatType rescaleMinMax(const FloatType& value, const FloatType& min, const FloatType& max);

FloatType interpolateLinear(const std::vector<double>& x, const std::vector<double>& y, const FloatType& xp);

void writeDataToCSV(const std::vector<FloatType>& data, const std::string& filename);
