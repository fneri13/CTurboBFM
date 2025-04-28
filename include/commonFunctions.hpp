#pragma once

#include "types.hpp"
#include "CFluid.hpp"

// compute the integral flux of an array2D of vectors through an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

// compute the total area of an array2D of surfaces
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces);

StateVector getEulerPrimitiveFromConservative(StateVector conservative);

StateVector getEulerConservativeFromPrimitive(StateVector primitive);

StateVector computeEulerFluxFromPrimitive(StateVector primitive, Vector3D surface, CFluid fluid);

StateVector computeEulerFluxFromConservative(StateVector conservative, Vector3D surface, CFluid fluid);
