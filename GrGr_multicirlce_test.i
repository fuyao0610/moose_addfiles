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
  op_num = 4
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 105
    grain_num = 4
    coloring_algorithm = bt
  [../]
[]

[ICs]
  [./PolycrystalICs]
#    [./PolycrystalColoringIC]
#      polycrystal_ic_uo = voronoi
#    [../]
     [./PolycrystalMultiCircleIC]  #syntax in PhaseFieldApp.C
      grain_num = 4
      circlespac = 400
      numtries = 1000
      radius = 150
      rand_seed = 98666
      radius_variation = 0.1
      radius_variation_type = uniform
      columnar_3D = false
     [../]
  [../]
#    [./PolycrystalRandomIC]
#      type = 0
#    [../]
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
  scheme = 'bdf2'

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 20
  nl_rel_tol = 1.0e-9
  start_time = 0.0
  num_steps = 100
  dt = 80.0
[]


[Outputs]
  file_base = multicircle
  exodus = true
[]
