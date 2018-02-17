/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POLYCRYSTALMULTICIRCLEICACTION_H
#define POLYCRYSTALMULTICIRCLEICACTION_H

#include "InputParameters.h"
#include "Action.h"

/**
 * Random multicircle polycrystal action
 */
class PolycrystalMultiCircleICAction : public Action
{
public:
  PolycrystalMultiCircleICAction(const InputParameters & params);

  virtual void act();

private:
  const unsigned int _op_num;
  const unsigned int _grain_num;

  const unsigned int _numtries;
  const unsigned int _rand_seed;

  const Real _circlespac;
  const Real _radius;
  const Real _radius_variation;
  //const MooseEnum _radius_variation_type;

  const std::string _var_name_base;

};

template <>
InputParameters validParams<PolycrystalMultiCircleICAction>();

#endif // POLYCRYSTALMULTICIRCLEICACTION_H
