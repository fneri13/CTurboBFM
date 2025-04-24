#pragma once

#include "types.hpp"
#include <iostream>
#include "CFluid.hpp"

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

// overloaded function, where only the summation is done
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces);

StateVector getEulerPrimitiveFromConservative(StateVector conservative);

StateVector getEulerConservativeFromPrimitive(StateVector primitive);

StateVector computeEulerFluxFromPrimitive(StateVector state, Vector3D surface, CFluid fluid);