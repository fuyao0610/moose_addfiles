/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ACGrGrPoly_solution.h"

template <>
InputParameters
validParams<ACGrGrPoly_solution>()
{
  InputParameters params = validParams<ACGrGrBase>();
  params.addClassDescription("Grain-Boundary model poly-crystaline interface Allen-Cahn Kernel");
  params.addRequiredParam<unsigned int>("op", "current order parameter");
  return params;
}

ACGrGrPoly_solution::ACGrGrPoly_solution(const InputParameters & parameters)
  : ACGrGrBase(parameters), 
    _gamma(getMaterialProperty<Real>("gamma_asymm")),
    _deltag(getMaterialProperty<Real>("deltag")),
    _op(getParam<unsigned int>("op"))
{
}

Real
ACGrGrPoly_solution::computeDFDOP(PFFunctionType type)
{
  // Sum all other order parameters
  Real SumEtaj = 0.0;
  Real SumEtaj_fenergy = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i) {
    if (i != _op)  {
      SumEtaj += (*_vals[i])[_qp] * (*_vals[i])[_qp];
      //      if (_op == 0) SumEtaj_fenergy += (*_vals[i])[_qp] * _deltag[_qp];  //liquid is represented by op = 0
      //if (i == 0) SumEtaj_fenergy -= (*_vals[i])[_qp] * _deltag[_qp];  
    }
  }

  // Calculate either the residual or Jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
    {
      const Real tgrad_correction =
          _grad_T ? _tgrad_corr_mult[_qp] * _grad_u[_qp] * (*_grad_T)[_qp] : 0.0;
      return _mu[_qp] *
	(_u[_qp] * _u[_qp] * _u[_qp] - _u[_qp] + 2.0 * _gamma[_qp] * _u[_qp] * SumEtaj) + SumEtaj_fenergy * _u[_qp] + tgrad_correction;
    }

    case Jacobian:
    {
      const Real tgrad_correction =
          _grad_T ? _tgrad_corr_mult[_qp] * _grad_phi[_j][_qp] * (*_grad_T)[_qp] : 0.0;
      return _mu[_qp] *
                 (_phi[_j][_qp] * (3.0 * _u[_qp] * _u[_qp] - 1.0 + 2.0 * _gamma[_qp] * SumEtaj)) + _phi[_j][_qp] * SumEtaj_fenergy  + tgrad_correction;
    }

    default:
      mooseError("Invalid type passed in");
  }
}

Real
ACGrGrPoly_solution::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i])
    {
      // Derivative of SumEtaj
      const Real dSumEtaj = 2.0 * (*_vals[i])[_qp] * _phi[_j][_qp];
      const Real dDFDOP = _mu[_qp] * 2.0 * _gamma[_qp] * _u[_qp] * dSumEtaj;

      return _L[_qp] * _test[_i][_qp] * dDFDOP;
    }

  return 0.0;
}
