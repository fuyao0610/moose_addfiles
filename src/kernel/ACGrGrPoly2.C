/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ACGrGrPoly2.h"
#include "Assembly.h"

template <>
InputParameters
validParams<ACGrGrPoly2>()
{
  InputParameters params = validParams<ACGrGrBase>();
  params.addClassDescription("Grain-Boundary model poly-crystaline interface Allen-Cahn Kernel");
  params.addRequiredParam<unsigned int>("op", "current order parameter");
  return params;
}

ACGrGrPoly2::ACGrGrPoly2(const InputParameters & parameters)
  : ACGrGrBase(parameters), 
    _gamma(getMaterialProperty<Real>("gamma_asymm")),
    //_sigma(getMaterialProperty<Real>("sigma")),
    //_l_GB(getMaterialProperty<Real>("l_GB")),
    _op(getParam<unsigned int>("op"))
{
}

Real
ACGrGrPoly2::computeDFDOP(PFFunctionType type)
{
  // Sum all other order parameters
  Real phi_j, lap_phi_j;  

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);

  Real SumEtaj = 0.0;

  for (unsigned int i = 0; i < _op_num; ++i) {
    if (_op != i) {

      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);

      //SumEtaj += (*_vals[i])[_qp] * (*_vals[i])[_qp];
      SumEtaj += phi_j * phi_j;
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
                 (phi_i * phi_i * phi_i - phi_i + 2.0 * _gamma[_qp] * phi_i * SumEtaj) +
             tgrad_correction;
    }

    case Jacobian:
    {
      const Real tgrad_correction =
          _grad_T ? _tgrad_corr_mult[_qp] * _grad_phi[_j][_qp] * (*_grad_T)[_qp] : 0.0;
      return _mu[_qp] *
                 (_phi[_j][_qp] * (3.0 * phi_i * phi_i - 1.0 + 2.0 * _gamma[_qp] * SumEtaj)) +
             tgrad_correction;
    }

    default:
      mooseError("Invalid type passed in");
  }
}

Real
ACGrGrPoly2::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real phi_j, lap_phi_j;  

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);

  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i] && _op != i) {
      // Derivative of SumEtaj
      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);
      const Real dSumEtaj = 2.0 * phi_j * _phi[_j][_qp];
      const Real dDFDOP = _mu[_qp] * 2.0 * _gamma[_qp] * phi_i * dSumEtaj;

      return _L[_qp] * _test[_i][_qp] * dDFDOP;
    }

  return 0.0;
}
