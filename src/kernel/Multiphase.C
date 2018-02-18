/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "Multiphase.h"
#include "Assembly.h"

template <>
InputParameters
validParams<Multiphase>()
{
  InputParameters params = validParams<ACGrGrBase>();
  params.addClassDescription("Grain-Boundary model poly-crystaline interface Allen-Cahn Kernel");
  params.addRequiredParam<unsigned int>("op", "current order parameter");
  return params;
}

Multiphase::Multiphase(const InputParameters & parameters)
  : ACGrGrBase(parameters), 
    _gamma(getMaterialProperty<Real>("gamma_asymm")),
    _kappa(getMaterialProperty<Real>("kappa_op")),
    _sigma(getMaterialProperty<Real>("sigma")),
    _l_GB(getMaterialProperty<Real>("l_GB")),
    _mu(getMaterialProperty<Real>("mu")),
    _deltag(getMaterialProperty<Real>("deltag")),
    _eps(getMaterialProperty<Real>("eps")),
    _W(getMaterialProperty<Real>("W")),
    _op(getParam<unsigned int>("op"))
{
}

Real
Multiphase::computeDFDOP(PFFunctionType type)
{

  Real Sumfenergy1 = 0.0;
  Real Sumfenergy2 = 0.0;

  Real dSumfenergy1 = 0.0;
  Real dSumfenergy2 = 0.0;

  Real phi_j;  

  Real SumPhi = 0.0;
  Real SumPhi_dg = 0.0;
  Real SumPhi_sq = 0.0;

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);

  // Sum all other order parameters
  for (unsigned int i = 0; i < _op_num; ++i) {
    if (_op != i) {
      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);

      SumPhi += phi_j;
      if (_op == 0) 
	SumPhi_dg += phi_j;
      else if (i == 0)
	SumPhi_dg -= phi_j;
      SumPhi_sq += phi_j * phi_j;
    }
  }


  // Calculate either the residual or Jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
    {
      Sumfenergy1 = 2 * _W[_qp] * (phi_i * SumPhi_sq - phi_i * phi_i * SumPhi);
      Sumfenergy2 = -6 * _deltag[_qp] * phi_i * SumPhi_dg;

      return Sumfenergy1 + Sumfenergy2;
    }

    case Jacobian:
    {
      dSumfenergy1 =  2 * _W[_qp] * (SumPhi_sq - 2 * phi_i * SumPhi) * _phi[_j][_qp];
      dSumfenergy2 =  -6 * _deltag[_qp] * _phi[_j][_qp] * SumPhi_dg;

      return dSumfenergy1 + dSumfenergy2; 
    }

    default:
      mooseError("Invalid type passed in");
  }
}

Real
Multiphase::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real dSumfenergy1 = 0.0;
  Real dSumfenergy2 = 0.0;
  Real phi_j;  

  Real SumPhi = 0.0;
  Real SumPhi_dg = 0.0;
  Real SumPhi_sq = 0.0;

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);

  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i] && _op != i) {

      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);

      dSumfenergy1 = 2 * _W[_qp] * (2 * phi_i * phi_j - phi_i * phi_i) * _phi[_j][_qp];
      if (_op > i)
	dSumfenergy2 = -6 * phi_i * _phi[_j][_qp] * _deltag[_qp];
      else
	dSumfenergy2 = 6 * phi_i * _phi[_j][_qp] * _deltag[_qp];

      const Real dDFDOP = dSumfenergy1 + dSumfenergy2;

      return _L[_qp] * _test[_i][_qp] * dDFDOP;
    }

  return 0.0;
}
