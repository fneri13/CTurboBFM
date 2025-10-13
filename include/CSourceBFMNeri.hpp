#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Thollet model source terms computation
class CSourceBFMNeri : public CSourceBFMBase {
public:

    CSourceBFMNeri(const Config &config, const CFluidBase &fluid, const CMesh &mesh, FlowSolution &conservativeVars, 
                    std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance);

    virtual ~CSourceBFMNeri() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) override;

    StateVector computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce);
    
    StateVector computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce);

    void computeStreamwiseCoefficient(size_t i, size_t j, size_t k);

    FloatType computeCompressibilityCorrection(const Vector3D& relativeVelocityCylindric, const StateVector& primitive);


private:

    FlowSolution &_conservativeSolution;
    
    FloatType _tangentialForce{0.0};
    
    CInputTable _inputTable;
    
    std::map<TurboPerformance, std::vector<FloatType>> &_turboPerformance;
    
    FloatType _scalingTurning{1.0}, _scalingLoss{1.0};
    
    size_t _trailingEdgeIndex{0}, _leadingEdgeIndex{0};
    
    Vector3D _viscousForceCylindrical{0, 0, 0};

    FloatType _streamwiseCoeff{0.0}, _mLeadingEdge{0.0}, _mTrailingEdge{0.0};
};
