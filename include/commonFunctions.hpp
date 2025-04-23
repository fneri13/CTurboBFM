#pragma once

#include "types.hpp"
#include <iostream>

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors);

// overloaded function, where only the summation is done
FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces);