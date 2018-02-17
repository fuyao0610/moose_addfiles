/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ACGRGRPOLY2_H
#define ACGRGRPOLY2_H

#include "ACGrGrBase.h"

//Forward Declarations
class ACGrGrPoly2;

template<>
InputParameters validParams<ACGrGrPoly2>();

/**
 * This kernel calculates the residual for grain growth for a single phase,
 * poly-crystal system. A single material property gamma_asymm is used for
 * the prefactor of the cross-terms between order parameters.
 */
class ACGrGrPoly2 : public ACGrGrBase
{
public:
  ACGrGrPoly2(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _gamma;
  const unsigned int _op;
};

#endif //ACGRGRPOLY2_H
