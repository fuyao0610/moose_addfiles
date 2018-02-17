/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef KKSMULTIMATERIALSACTION_H
#define KKSMULTIMATERIALSACTION_H

#include "InputParameters.h"
#include "Action.h"

/**
 * Automatically generates all variables, Kernels, and Materials to ensure the
 * correct derivatives of the elastic free energy in a non-split Cahn-Hilliard
 * simulation are assembled.
 */
class KKSMultiMaterialsAction : public Action
{
public:
  KKSMultiMaterialsAction(const InputParameters & params);

  virtual void act();

private:
  unsigned int _op_num;
  std::vector<Real> _eq_c;
};

template <>
InputParameters validParams<KKSMultiMaterialsAction>();

#endif // KKSMULTIMATERIALSACTION_H
