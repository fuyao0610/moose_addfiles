/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef KKSVARIABLESACTION_H
#define KKSVARIABLESACTION_H

#include "InputParameters.h"
#include "Action.h"

/**
 * Automatically generates all variables to model a polycrystal with op_num orderparameters
 */
class KKSVariablesAction : public Action
{
public:
  KKSVariablesAction(const InputParameters & params);

  virtual void act();

private:
  const unsigned int _op_num;
  //const std::string _var_name_base;
};

template <>
InputParameters validParams<KKSVariablesAction>();

#endif // KKSVARIABLESACTION_H
