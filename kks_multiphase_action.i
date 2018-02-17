#
# This test is for the 3-phase KKS model
#

[GlobalParams]
  op_num = 4
  var_name_base = eta
  xmin = 0
  xmax = 40
  ymin = 0
  ymax = 40
  zmin = 0
  zmax = 0
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  nz = 0
  elem_type = QUAD4
[]


[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Variables]
  # concentration
  [./KKSVariables]
  [../]
[]

[AuxVariables]
   [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Materials]
  # simple toy free energies
  [./KKSMaterials]
    eq_concentration = '0.2 0.6' 
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L   kappa  D'
    prop_values = '0.7 1.0    1'
  [../]
[]

[Kernels]
  [./KKSMultiphase]
  [../]
[]

[ICs]
  [./KKSMultiCircleIC]
     circlespac = 14
     int_width = 2
     rand_seed = 1772
     avg_radius = 5
     radius_variation = 0.05
     twophase_eq_concentration = '0.45 0.6'	
  [../]	
[]


[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'F1 F2 F3'
    hj_names = 'h1 h2 h3'
    gj_names = 'g1 g2 g3'
    variable = Energy
    w = 1
    interfacial_vars =  'eta1  eta2  eta3'
    kappa_names =       'kappa kappa kappa'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-11

  start_time = 0.0
  end_time = 400

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]
[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    full = true
  [../]
  [./mydebug]
    type = FDP
    full = true
  [../]
[]

[Outputs]
  exodus = true
[]
