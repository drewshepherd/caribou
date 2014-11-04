###################################################
##2D with Q_fission heat generation, contact and burnup-dependent heat conduction
##taken from heat_source_bar.i
###################################################

[GlobalParams]
	family = LAGRANGE
	order = SECOND
	#disp_x = disp_x
	#disp_y = disp_y
[]


[Problem]
	coord_type = RZ
[]


[Mesh]
	file = rmctests/2D/2D_simple_v1.e
	displacements = 'disp_x disp_y'
	patch_size = 1000
[]


[Functions]
	[./temperature]
		type = ConstantFunction
		value = 400
	[../]
	[./linear_power_function]
		type = PiecewiseLinear
		x = '0 10'
		y = '0 43000'
	[../]
[]


[AuxVariables]
	[./avg_burnup]
		block = pellets
	[../]
	[./burnup]
		block = pellets
		order = FIRST ##################################why is this a constant and monomial??###################3
		family = MONOMIAL
	[../]
[]


[Variables]
	[./temp]
		initial_condition = 400.0
	[../]
	[./disp_x]
  [../]
  [./disp_y]
  [../]
[]


[SolidMechanics]
  [./solid]
  	disp_r = disp_x
  	disp_z = disp_y
	  temp = temp
  [../]
[]

[AuxKernels]
	[./avg_burnup_aux]
		type = AverageBurnupAux
		linear_power = linear_power_function
		pellet_length = 15e-3
		initial_fuel_density = 1e9
		number_pellets = 1.5
		pellet_volume = 1.798368e-6
		ratio = .925
		variable = avg_burnup
	[../]
	[./burnup_aux]
		type = BurnupAux
		ratio = .925
		variable = burnup
		execute_on = timestep
	[../]
[]

[Kernels]
	[./heat]
		type = HeatConduction
		variable = temp
	[../]
	[./heat_dt]
		type = HeatConductionTimeDerivative		
		variable = temp
	[../]
	[./heatsource]
		type = FissionHeatKernel
		variable = temp
		block = pellets
	[../]
[]


[BCs]
	[./sheath_out]
		type = ConvectiveFluxFunction
    boundary = sheath_out
    variable = temp
		T_infinity = temperature 
		coefficient = 5e4
	[../]

	[./pelletsz2]
    type = DirichletBC
    variable = disp_y
    boundary = top_pellet2
    value = 0
  [../]
#	[./pelletsz]
#    type = DirichletBC
#    variable = disp_y
#    boundary = bot_land_point
#    value = 0
#  [../]

  [./pelletsr]
    type = DirichletBC
    variable = disp_x
    boundary = pellets_in
    value = 0
  [../]
	
  [./sheath_z]
    type = DirichletBC
    variable = disp_y
    boundary = 'sheath_top'
    value = 0
  [../]	

	[./Pressure]
		[./coolant_pressure]
			disp_x = disp_x
			disp_y = disp_y
			boundary = 'sheath_out'
			factor = 0
  	[../]
	[../]
[]


[Contact]
  [./pellet_to_sheath]
		slave = pellets_out
    master = sheath_in
		formulation = penalty
		model = frictionless
		penalty = 1e8
		disp_x = disp_x
    disp_y = disp_y
  [../]
[]


[ThermalContact]
  [./pellet_to_sheath]
		type = GapHeatTransfer
#		gap_conductivity= 0.15
    slave = pellets_out
    master = sheath_in
    variable = temp
		quadrature = true
  [../]
[]


[Materials]
	 [./hcm]
    type = HeatConductionMaterial_burnup_temp
    block = pellets
    specific_heat = 320
    thermal_conductivity = 3
		temp = temp
		burnup = burnup
  [../]
    [./hcm_sheath]
    type = HeatConductionMaterial
    block = sheath
    specific_heat = 330
    thermal_conductivity = 16
		temp = temp
  [../]
  [./sm_materials_pellet]
    type = Elastic
		block = 'pellets'
		disp_r = disp_x
		disp_z = disp_y
		temp = temp
		thermal_expansion = 1e-5
		youngs_modulus = 100e9
		poissons_ratio = 0.3
		stress_free_temperature = 300
#		formulation = AxisymmetricRZ
  [../]
  [./sm_materials_sheath]
    type = Elastic
		block = 'sheath'
		disp_r = disp_x
		disp_z = disp_y
		temp = temp
		thermal_expansion = 0
		youngs_modulus = 10e9
		poissons_ratio = 0.3
		stress_free_temperature = 400
#		formulation = AxisymmetricRZ
  [../]
	[./pellet_density]
    type = Density
		block = pellets
		disp_r = disp_x
		disp_z = disp_y
		density = 1.065e4
  [../]
	[./sheath_density]
    type = Density
		block = sheath
		disp_r = disp_x
		disp_z = disp_y
		density = 6.5e3
  [../]
	[./FissionHeatMaterial]
		type = FissionHeatMaterial
		block = pellets
		initial_fuel_density = 10e3
		initial_qfission = 0
		enrichment = 5
		burnup_avg = avg_burnup
		pellet_radius = 6e-3
		linear_power = linear_power_function
	[../]
[]


[Dampers]
	[./MaxIncrement]
		type = MaxIncrement
		variable = temp
		max_increment = 50
	[../]	
	[./MaxIncrementDispX]
		type = MaxIncrement
		variable = disp_x
		max_increment = 1e-5
	[../]
	[./MaxIncrementDispY]
		variable = disp_y
		max_increment = 1e-5
		type = MaxIncrement
	[../]
[]


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
  end_time = 1e2
	num_steps = 5000

	[./TimeStepper]
		type = IterationAdaptiveDT
		dt = 1
		cutback_factor = 0.5
		growth_factor = 2
		optimal_iterations = 10
		iteration_window = 2 #0.2*optimal_iterations
		linear_iteration_ratio = 100
	[../]
[]


[Outputs]
  file_base = rmctests/2D/2D_simple_4_burnupaux_out
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
#   linear_residuals = true
  [../]
[]


#[Debug]
#	show_var_residual_norms = true
#[]
