/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef LATENTHEAT_H
#define LATENTHEAT_H

#include "Kernel.h"

// Forward Declaration
class Latentheat;

template <>
InputParameters validParams<Latentheat>();

/**
 * Kernel which implements the convective term in the transient heat
 * conduction equation, and provides coupling with the Darcy pressure
 * equation.
 */
class Latentheat : public Kernel
{
public:
  Latentheat(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

#  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// The gradient of pressure
  const Real  & _phase_dot;

  /// Coupling identifier for the pressure.  This is used to uniquely
  /// identify a coupled variable
  unsigned int _phase_var;

  /// These references will be set by the initialization list so that
  /// values can be pulled from the Material system.
  const MaterialProperty<Real> & _latentheat;
};

#endif // LATENTHEAT_H
