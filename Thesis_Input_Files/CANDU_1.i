########################################################
##CANDU geometry, 1/4 rod in the first ring
########################################################

[GlobalParams]
  family = LAGRANGE
  order = SECOND
[]

[Problem]
	coord_type = XYZ
[]

##Trelis geometry
[Mesh]
  file = Thesis_Geometry/CANDU_cross_section.e
  displacements = 'disp_x disp_y disp_z'
  patch_size = 1000
[] ##Mesh

##Functions to be used in the Kernels, AuxKernels and Materials
[Functions]
  [./linear_power_function_time] ##scale this for volume and then increase with increasing radius, find a way to make this bundle power the same as a CANDU,
    type = PiecewiseLinear
    x = '0 10'
    y = '0 43000' ##[W/m]
  [../]
	[./linear_power_function_burnup]
		type = PiecewiseLinear
		x = '0'
		y = '0'
	[../]
  [./pressure_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1'
  [../]
[] ##Functions

##Variables to be solved for in the AuxKernels
[AuxVariables]
  [./avg_burnup] ##solves for the average burnup - only used in FissionHeatMaterial
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./burnup] ##solves for the local burnup - used everywhere except FissionHeatMaterial
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./burnup_dt] ##solves for the local burnup_dt
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./vstrain] ##solves for the volumetric strain
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./k_pellets]
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./k_pellets_dT]
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
  [./k_sheath]
    order = FIRST
    family = MONOMIAL
    block = sheath
  [../]
  [./k_sheath_dT]
    order = FIRST
    family = MONOMIAL
    block = sheath
  [../]
  [./test]
    order = FIRST
    family = MONOMIAL
    block = pellet
  [../]
[] ##AuxVariables

##Variables to be solved for in the differential equations (ODE/PDE) in the Kernels
[Variables]
  [./temp] ##temperature variable
    initial_condition = 550
    block = 'pellet sheath'
  [../]
  [./disp_x] ##x-displacement variable
  [../]
  [./disp_y] ##y-displacement variable
  [../]
  [./disp_z] ##z-displacement variable
  [../]
  [./densification_fraction] ##densification fraction variable - a percentage less than 60%
    initial_condition = 0
    block = pellet
  [../]
[] ##Variables

##Solve for the AuxVariables - these do not solve ODE/PDE
[AuxKernels]
  [./avg_burnup_aux] ##solves for avg_burnup auxvariable
    type = AverageBurnupAux
    variable = avg_burnup
    initial_density = 1.065e4
    execute_on = timestep_begin
    model_wrt_time = true
    linear_power_time = linear_power_function_time
    model_wrt_burnup = false
    linear_power_burnup = linear_power_function_burnup
  [../]
  [./burnup_aux] ##solves for burnup auxvariable
    type = BurnupAux
    variable = burnup
    execute_on = timestep_begin
    model_burnup = true
  [../]
  [./burnup_dt_aux] ##solves for burnup auxvariable
    type = Burnup_dtAux
    variable = burnup_dt
    execute_on = timestep_begin
    model_burnup_dt = true
  [../]
  [./vstrain_aux] ##solves for burnup auxvariable
    type = VStrainAux
    variable = vstrain
    densification_fraction = densification_fraction
    execute_on = timestep_begin
    initial_porosity = 0.0286
    model_vstrain = true
  [../]
  [./k_pellets_aux] ##solves for burnup auxvariable
    type = ThermalConductivityPelletsAux
    variable = k_pellets
    temp = temp
    burnup = burnup
    execute_on = timestep_begin
    model_k = true
  [../]
  [./k_pellets_dT_aux] ##solves for burnup auxvariable
    type = ThermalConductivity_dTPelletsAux
    variable = k_pellets_dT
    temp = temp
    burnup = burnup
    execute_on = timestep_begin
    model_k_dT = true
  [../]
  [./k_sheath_aux] ##solves for burnup auxvariable
    type = ThermalConductivitySheathAux
    variable = k_sheath
    temp = temp
    execute_on = timestep_begin
    model_k = true
  [../]
  [./k_sheath_dT_aux] ##solves for burnup auxvariable
    type = ThermalConductivity_dTSheathAux
    variable = k_sheath_dT
    temp = temp
    execute_on = timestep_begin
    model_k_dT = true
  [../]
  [./test_aux]
    type = TestAux
    variable = test
    execute_on = timestep_begin
    model_density = true
    model_q_fission = false
    model_q_fission_old = false
    model_ratio = false
    model_pellet_radius = false
    model_porosity = false
    model_enrichment = false
    model_SFP = false
    model_GFP = false
  [../]
[] ##AuxKernels

##Solve for the Variables
[Kernels]
  [./HC] ##solves for the thermal conductivity - couples to Heat_dt, HeatSource and Densification_Fraction
    type = HeatConduction
    variable = temp
  [../]
  [./Heat_dt] ##solves for the time rate of change of the thermal conductivity - couples to HC, HeatSource and Densification_Fraction
    type = HeatConductionTimeDerivative
    variable = temp
  [../]
  [./HeatSource] ##solves for the temperature of the pellet given the heat produced from nuclear fission (uses the FissionHeatMaterial) - couples to HC, Heat_dt and Densification_Fraction
    type = FissionHeatKernel
    variable = temp
    block = pellet
  [../]
  [./Densification_Fraction_RHS] ##solves for F_porosity - couples to Denisification_Fraction_LHS
    type = DensificationFractionKernel
    variable = densification_fraction ##this is a decimal value
    temp = temp
    burnup_dt = burnup_dt
    model_densification_fraction = true
  [../]
  [./Densification_Fraction_LHS] ##solves for F_porosity - couples to Denisification_Fraction_RHS
    type = TimeDerivative
    variable = densification_fraction ##this is a decimal value
  [../]
[] ##Kernels

##Allows for solid mechanics to be used
[SolidMechanics]
  [./solid]
    temp = temp
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[] ##SolidMechanics

[BCs]
  [./Sheath_out] ##Specifies the coolant temperature and the rate of heat flux out of the bundle
    type = ConvectiveFluxFunction
    variable = temp
    boundary = 'sheath_out'
    T_infinity = 550
    coefficient = 5e4
  [../]
  [./Y] ##Constrains the midline_x of the pellet and sheath to not move vertically
    type = DirichletBC
    variable = disp_y
    boundary = 'pellet_x sheath_x'
    value = 0
  [../]
  [./X] ##Constrains the midline_y of the pellet and sheath to not move horizontally
    type = DirichletBC
    variable = disp_x
    boundary = 'pellet_y sheath_y'
    value = 0
  [../]
  [./Z] ##Constrains the pellet and sheath to not move vertically
    type = DirichletBC
    variable = disp_z
    boundary = 'pellet_bot sheath_bot' ##allows the pellet to expand in one direction in Z
    value = 0
  [../]
  [./Pressure] ##applies the 10MPa pressure to the sheath
    [./coolant_pressure]
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      boundary = sheath_out
      factor = 1e7
      function = pressure_ramp
    [../]
  [../]
[] ##BCs

##Specifies the properties when the pellet-sheath contact each other - specifies the penalty factor coefficient
[Contact]
  [./pellet_to_sheath_contact]
    slave = pellet_out
    master = sheath_in
    formulation = penalty
    model = frictionless
    penalty = 1e8
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[] ##Contact

##Allows heat to transfer over the gap between the pellet-sheath
[ThermalContact]
  [./pellet_to_sheath]
    type = GapHeatTransfer
    variable = temp
    slave = pellet_out
    master = sheath_in
    quadrature = true
    gap_conductivity= 0.15
  [../]
[] ##ThermalContact

[Materials]
  [./HCM_Pellet] ##Determines the thermal properties of the pellet - thermal conductivity, specific heat, porosity. If these models are set to false, then default values of 1 are given
    type = PelletThermalMaterial
    temp = temp
    burnup = burnup
    burnup_dt = burnup_dt
    densification_fraction = densification_fraction
    k_pellets = k_pellets
    k_pellets_dT = k_pellets_dT
    block = pellet
    initial_porosity = 0.0286
    model_thermal_conductivity = true
    model_specific_heat = true
    model_porosity = true ##this causes no change in density, as it is not yet coupled to density
    model_alpha = true
    model_SFP = true
    model_GFP = true
    display_values = false
  [../]
  [./HCM_Sheath] ##Determines the thermal properties of the sheath - thermal conductivity, specific heat, thermal expansion
    type = SheathThermalMaterial
    temp = temp
    block = sheath
    k_sheath = k_sheath
    k_sheath_dT = k_sheath_dT
    model_thermal_conductivity = true
    model_specific_heat = true
  [../]
  [./Mechanical_Pellets] ##Determines the solid mechanical properties of the pellet - stresses, strains
    type = StrainMaterial
    temp = temp
    burnup = burnup
    burnup_dt = burnup_dt
    vstrain = vstrain
    block = pellet
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    poissons_ratio = 0.3
    youngs_modulus = 1.8e11
    model_thermal_expansion = true
    model_youngs_modulus = true
    display_values = false
  [../]
  [./Mechanical_Sheath] ##Determines the solid mechanical properties of the sheath - stresses, strains
    type = SheathMechanicalMaterial
    temp = temp
    block = sheath
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    poissons_ratio = 0.3
    youngs_modulus = 1.8e11
    relative_tolerance = 1e-4
    absolute_tolerance = 1e-9
    max_its = 14
    output_iteration_info = false
    model_diffusional_creep = false
    model_thermal_expansion = true
    model_youngs_modulus = true
    display_values = false
  [../]
  [./Pellet_density ] ##Specifies the density of the pellets
    type = Density
    block = pellet
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    density = 1.065e4
  [../]
  [./Sheath_density] ##Specifies the density of the sheath
    type = Density
    block = sheath
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    density = 6.5e3
  [../]
  [./FissionHeat] ##Determines the amount of heat produced by nuclear fission
    type = FissionHeatMaterial
    burnup_avg = avg_burnup
    block = pellet
    linear_power = linear_power_function_time
    initial_fuel_density = 1.065e4
    initial_qfission = 0
    enrichment = 1
    pellet_radius = 6.118e-3
    ratio = 0.925
		initial_fuel_area = 1.1767e-4 ##in m^2, for a full rod
    is_3D = true
    model_Qfission = true
		model_plate_fuel = false
  [../]
[] ##Materials


[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      8'

  line_search = 'none'

  #petsc options from previous example
  #	petsc_options_iname = '-pc_type -pc_hypre_type'
  #  petsc_options_value = 'hypre boomeramg'

  # controls for linear iterations
  l_max_its = 100
  l_tol = 1e-3 #3

  # controls for nonlinear iterations
  nl_max_its = 15
  nl_abs_tol = 5e-9 #10
  nl_rel_tol = 1e-5 #5

  #time control
  start_time = 0.0
  #  dt = 10.0
  end_time = 2.5056e7 #290 days in the reactor
  num_steps = 5000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = .25
    cutback_factor = 0.5
    growth_factor = 2
    optimal_iterations = 10
    iteration_window = 2 #0.2*optimal_iterations
    linear_iteration_ratio = 100
  [../]
[] ##Executioner


[Outputs]
  file_base = Thesis_Paraview/CANDU_1_out
  output_initial = true
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[] ##Outputs


[Debug]
  show_var_residual_norms = false
[] ##Debug
