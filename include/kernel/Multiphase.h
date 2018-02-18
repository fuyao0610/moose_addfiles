/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "ACGrGrBase.h"

// Forward Declarations
class Multiphase;

template <>
InputParameters validParams<Multiphase>();

/**
 * This kernel calculates the residual for grain growth for a single phase,
 * poly-crystal system. A single material property gamma_asymm is used for
 * the prefactor of the cross-terms between order parameters.
 */
class Multiphase : public ACGrGrBase
{
public:
  Multiphase(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  const MaterialProperty<Real> & _gamma;  
  const MaterialProperty<Real> & _kappa;  
  const MaterialProperty<Real> & _sigma;
  const MaterialProperty<Real> & _l_GB;
  const MaterialProperty<Real> & _mu;  
  const MaterialProperty<Real> & _deltag;
  const MaterialProperty<Real> & _eps;  
  const MaterialProperty<Real> & _W;
  const unsigned int _op;

};

#endif // MULTIPHASE_H
