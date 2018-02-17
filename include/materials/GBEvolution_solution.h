/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GBEVOLUTION_SOLUTION_H
#define GBEVOLUTION_SOLUTION_H

#include "GBEvolutionBase.h"

// Forward Declarations
class GBEvolution_solution;

template <>
InputParameters validParams<GBEvolution_solution>();

/**
 * Grain boundary energy parameters for isotropic uniform grain boundary energies
 */
class GBEvolution_solution : public GBEvolutionBase
{
public:
  GBEvolution_solution(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  Real _GBEnergy; //grain boundary surface energy 
  Real _Tm;  //melting temperature
  Real _LatentHeat;  //J/m^3

  MaterialProperty<Real> &  _deltag;
  MaterialProperty<Real> &  _eps;
  MaterialProperty<Real> &  _mobility;
  MaterialProperty<Real> &  _W;
};

#endif // GBEVOLUTION_SOLUTION_H
