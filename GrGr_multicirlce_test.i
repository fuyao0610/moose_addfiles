[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  nz = 0
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 5
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalMultiCircleIC]  #syntax in PhaseFieldApp.C
      grain_num = 5
      circlespac = 400
      numtries = 1000
      radius = 150
      rand_seed = 28872
      radius_variation = 0.1
      radius_variation_type = uniform
      columnar_3D = false
    [../]
  [../]
#    [./PolycrystalRandomIC]
#      type = 0
#    [../]
  [../]
[]


[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./PolycrystalMultiphaseKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution_solution
    T = 1000 # K
    wGB = 60  # nm
    GBMobility = 1.2e-8 #m^4/(Js) from Schoenfelder 1997
    # Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    LatentHeat = 2.31E9  #in J/m^3
    Tm = 1710 #K
  [../]
[]

[Postprocessors]
  active = ''
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
[]

[Preconditioning]
  active = ''
  [./SMP]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient
  scheme = 'implicit-euler'

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  start_time = 0.0
  num_steps = 5000
  dt = 0.02
[]


[Outputs]
  file_base = multicircle
  exodus = true
[]
