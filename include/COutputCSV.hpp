#pragma once

#include "COutputBase.hpp"

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class COutputCSV : public COutputBase {
    
    public:
        COutputCSV(const Config &config, const CMesh &mesh, const FlowSolution &solution) : COutputBase(config, mesh, solution) {};
        
        virtual ~COutputCSV() = default;

        virtual void writeSolution() override;

};