/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "GBEvolution_solution.h"

template <>
InputParameters
validParams<GBEvolution_solution>()
{
  InputParameters params = validParams<GBEvolutionBase>();
  params.addRequiredParam<Real>("GBenergy", "Grain boundary energy in J/m^2");
  params.addRequiredParam<Real>("LatentHeat", "Latent heat in J/m^3");
  params.addRequiredParam<Real>("Tm", "Melting temperature in K");
  return params;
}

GBEvolution_solution::GBEvolution_solution(const InputParameters & parameters)
  : GBEvolutionBase(parameters), 
    _GBEnergy(getParam<Real>("GBenergy")),
    _Tm(getParam<Real>("Tm")),
    _LatentHeat(getParam<Real>("LatentHeat")), 
    _deltag(declareProperty<Real>("deltag")),
    _eps(declareProperty<Real>("eps")),
    _mobility(declareProperty<Real>("mobility")),
    _W(declareProperty<Real>("W"))
{
}

void
GBEvolution_solution::computeQpProperties()
{
  // eV/nm^2
  _sigma[_qp] = _GBEnergy * _JtoeV * (_length_scale * _length_scale);
  
  Real _latentheat_convert = _LatentHeat;
  _latentheat_convert *= _JtoeV * (_length_scale * _length_scale * _length_scale);
  _deltag[_qp] = _latentheat_convert * (_T[_qp] - _Tm) / _Tm; 
  _eps[_qp] = 6 * _sigma[_qp] * _l_GB[_qp];
  _W[_qp] = 3 * _sigma[_qp] / _l_GB[_qp];
  _mobility[_qp] = _M_GB[_qp] / _l_GB[_qp];;
 
  GBEvolutionBase::computeQpProperties();
}
