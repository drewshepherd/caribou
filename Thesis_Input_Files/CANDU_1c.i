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
  file = Thesis_Geometry/pellet.e
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
[] ##BCs


[Materials]
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
  file_base = Thesis_Paraview/CANDU_1c_out
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = false
  [../]
[] ##Outputs


[Debug]
  show_var_residual_norms = false
[] ##Debug
