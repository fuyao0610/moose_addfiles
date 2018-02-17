/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PolycrystalMultiCircleICAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"
#include "PolycrystalICTools.h"

template <>
InputParameters
validParams<PolycrystalMultiCircleICAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Random multicircle polycrystal action");
  params.addRequiredParam<unsigned int>("op_num", "number of order parameters to create");
  params.addRequiredParam<unsigned int>(
      "grain_num", "number of grains to create, if it is going to greater than op_num");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");

  params.addRequiredParam<Real>("circlespac",
                                "minimum spacing of circles, measured from center to center");
  params.addParam<unsigned int>("numtries", 1000, "The number of tries");
  params.addParam<unsigned int>("rand_seed", 12444, "The random seed");
  params.addRequiredParam<Real>("radius", "Mean radius value for the circles");
  params.addParam<Real>("radius_variation",
                        0.0,
                        "Plus or minus fraction of random variation in "
                        "the grain radius for uniform, standard "
                        "deviation for normal");
  MooseEnum rand_options("uniform none", "none");
  params.addParam<MooseEnum>("radius_variation_type",
                             rand_options,
                             "Type of distribution that random circle radii will follow");
  params.addParam<bool>(
      "columnar_3D", false, "3D microstructure will be columnar in the z-direction?");
  return params;
}

PolycrystalMultiCircleICAction::PolycrystalMultiCircleICAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _grain_num(getParam<unsigned int>("grain_num")),
    _numtries(getParam<unsigned int>("numtries")),
    _rand_seed(getParam<unsigned int>("rand_seed")),
    _circlespac(getParam<Real>("circlespac")),
    _radius(getParam<Real>("radius")),
    _radius_variation(getParam<Real>("radius_variation")),
    _var_name_base(getParam<std::string>("var_name_base"))
    //_radius_variation_type(getParam<MooseEnum>("radius_variation_type"))
{
}

void
PolycrystalMultiCircleICAction::act()
{
#ifdef DEBUG
  Moose::err << "Inside the PolycrystalMultiCircleICAction Object\n";
#endif

  // Loop through the number of order parameters
  for (unsigned int op = 0; op < _op_num; op++)
  {
    // Set parameters for BoundingBoxIC
    InputParameters poly_params = _factory.getValidParams("PolycrystalMultiCircleIC");
    poly_params.set<VariableName>("variable") = _var_name_base + Moose::stringify(op);
    poly_params.set<unsigned int>("op_num") = _op_num;
    poly_params.set<unsigned int>("grain_num") = _grain_num;
    poly_params.set<unsigned int>("op_index") = op;

    poly_params.set<Real>("circlespac") = _circlespac;
    poly_params.set<Real>("numtries") = _numtries;
    poly_params.set<Real>("radius") = _radius;
    poly_params.set<unsigned int>("rand_seed") = _rand_seed;
    poly_params.set<Real>("radius_variation") = _radius_variation;
    poly_params.set<MooseEnum>("radius_variation_type") = getParam<MooseEnum>("radius_variation_type");
    poly_params.set<bool>("columnar_3D") = getParam<bool>("columnar_3D");

    // Add initial condition
    _problem->addInitialCondition(
        "PolycrystalMultiCircleIC", "PolycrystalCircleIC_" + Moose::stringify(op), poly_params);
  }
}
