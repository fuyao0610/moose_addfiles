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
    _op(getParam<unsigned int>("op")),
    _second_phi(_assembly.secondPhi()),
    _laplac_vals(_op_num)
{

  for (unsigned int i = 0; i < _op_num; ++i)
      _laplac_vals[i] = &coupledSecond("v", i);
 
}

Real
Multiphase::computeDFDOP(PFFunctionType type)
{

  Real Sumintface = 0.0;
  Real Sumfenergy1 = 0.0;
  Real Sumfenergy2 = 0.0;

  Real dSumintface = 0.0;
  Real dSumfenergy1 = 0.0;
  Real dSumfenergy2 = 0.0;

  Real phi_j, lap_phi_j;  

  Real SumPhi = 0.0;
  Real SumPhi_dg = 0.0;
  Real SumPhi_sq = 0.0;
  Real SumPhi_laplac = 0.0;

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);
  Real lap_phi_i = (*_laplac_vals[_op])[_qp].tr();

  // Sum all other order parameters
  for (unsigned int i = 0; i < _op_num; ++i) {
    if (_op != i) {
      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);
      lap_phi_j = (*_laplac_vals[i])[_qp].tr();

      SumPhi += phi_j;
      if (_op > i) 
	SumPhi_dg += phi_j;
      else 
	SumPhi_dg -= phi_j;
      SumPhi_sq += phi_j * phi_j;
      SumPhi_laplac += lap_phi_j;
    }
  }


  // Calculate either the residual or Jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
    {
      Sumintface = _eps[_qp] * (phi_i * SumPhi_laplac - lap_phi_i * SumPhi);
      Sumfenergy1 = 2 * _W[_qp] * (phi_i * SumPhi_sq - phi_i * phi_i * SumPhi);
      Sumfenergy2 = -6 * _deltag[_qp] * phi_i * SumPhi_dg;

      return Sumintface + Sumfenergy1 + Sumfenergy2;
    }

    case Jacobian:
    {
      dSumintface = _eps[_qp] * (_phi[_j][_qp] * SumPhi_laplac - _second_phi[_j][_qp].tr() * SumPhi);
      dSumfenergy1 =  2 * _W[_qp] * (SumPhi_sq - 2 * phi_i * SumPhi) * _phi[_j][_qp];
      dSumfenergy2 =  -6 * _deltag[_qp] * _phi[_j][_qp] * SumPhi_dg;

      return dSumintface + dSumfenergy1 + dSumfenergy2; 
    }

    default:
      mooseError("Invalid type passed in");
  }
}

Real
Multiphase::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real dSumintface = 0.0;
  Real dSumfenergy1 = 0.0;
  Real dSumfenergy2 = 0.0;
  Real phi_j, lap_phi_j;  

  Real SumPhi = 0.0;
  Real SumPhi_dg = 0.0;
  Real SumPhi_sq = 0.0;
  Real SumPhi_laplac = 0.0;

  Real phi_i = fmax(fmin(_u[_qp], 1), 0);
  Real lap_phi_i = (*_laplac_vals[_op])[_qp].tr();

  for (unsigned int i = 0; i < _op_num; ++i)
    if (jvar == _vals_var[i] && _op != i) {

      phi_j = fmax(fmin((*_vals[i])[_qp], 1), 0);
      lap_phi_j = (*_laplac_vals[i])[_qp].tr();

      dSumintface = _eps[_qp] * (phi_i * _second_phi[_j][_qp].tr() - _phi[_j][_qp] * lap_phi_i);
      dSumfenergy1 = 2 * _W[_qp] * (2 * phi_i * phi_j - phi_i * phi_i) * _phi[_j][_qp];
      if (_op > i)
	dSumfenergy2 = -6 * phi_i * _phi[_j][_qp] * _deltag[_qp];
      else
	dSumfenergy2 = 6 * phi_i * _phi[_j][_qp] * _deltag[_qp];

      const Real dDFDOP = dSumintface + dSumfenergy1 + dSumfenergy2;

      return _L[_qp] * _test[_i][_qp] * dDFDOP;
    }

  return 0.0;
}
