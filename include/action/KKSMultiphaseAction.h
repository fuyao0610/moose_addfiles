/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef KKSMULTIPHASEACTION_H
#define KKSMULTIPHASEACTION_H

#include "Action.h"

/**
 * Action that sets up ACGrGrPoly, ACInterface, TimeDerivative, and ACGBPoly
 * kernels.
 */
class KKSMultiphaseAction : public Action
{
public:
  KKSMultiphaseAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  //const std::string _var_name_base;

};

template <>
InputParameters validParams<KKSMultiphaseAction>();

#endif // KKSMULTIPHASEACTION_H
