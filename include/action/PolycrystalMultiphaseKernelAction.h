/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POLYCRYSTALMULTIPHASEKERNELACTION_H
#define POLYCRYSTALMULTIPHASEKERNELACTION_H

#include "Action.h"

/**
 * Action that sets up Multiphase, MultiInterface, and TimeDerivative
 * kernels.
 */
class PolycrystalMultiphaseKernelAction : public Action
{
public:
 PolycrystalMultiphaseKernelAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _var_name_base;

  /// kernels are implicit?
  const bool _implicit;
};

template <>
InputParameters validParams<PolycrystalMultiphaseKernelAction>();

#endif // POLYCRYSTALMULTIPHASEKERNELACTION_H
