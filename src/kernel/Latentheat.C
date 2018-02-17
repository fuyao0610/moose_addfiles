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

#include "Latentheat.h"

template <>
InputParameters
validParams<Latentheat>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("phase_variable", "The variable representing the phase.");

  return params;
}

Latentheat::Latentheat(const InputParameters & parameters)
  : Kernel(parameters),

    // Couple to the gradient of the pressure
    _phase_dot(coupledDot("phase_variable")),

    // Grab necessary material properties
    _latentheat(getMaterialProperty<Real>("latent_heat")),
{
}

Real
Latentheat::computeQpResidual()
{
  return _latentheat[_qp] * _phase_dot[_qp] * _test[_i][_qp];
}

Real
Latentheat::computeQpJacobian()
{
  return 0.0;
}

/*Real
Latentheat::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _pressure_var)
  {
    RealVectorValue superficial_velocity =
        _porosity[_qp] * -(_permeability[_qp] / _viscosity[_qp]) * _grad_phi[_j][_qp];
    return _heat_capacity[_qp] * superficial_velocity * _grad_u[_qp] * _test[_i][_qp];
  }

  return 0.0;
  }*/


