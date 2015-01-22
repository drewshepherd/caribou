###################################################
##2D with heat genration that depends on time 
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
	[./heat_gen]
		type = ParsedFunction
		value = 20e10*(6.44e-3-x) #there is no heat generation at sheath_out, and maximum generation at pellets_in
		#value = 3e8
	[../]
	[./body]
		type = ParsedFunction
		value = 1e1*x
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
	[./burnup]
		initial_condition = 1
	[../]
[]


[SolidMechanics]
  [./solid]
  	disp_r = disp_x
  	disp_z = disp_y
	  temp = temp
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
		type = HeatSource
		block = pellets
		variable = temp
		function = heat_gen	
	[../]
	[./burnupsource] ## [MWh/kgU], fraction of U atoms that have been fissioned
		type = BodyForce
		variable = burnup
		block = 'pellets sheath'
		value = 1		
		function = body
		use_displaced_mesh = 0
	[../]
	[./burnup_dt]
		type = TimeDerivative
		block = 'pellets sheath'
		variable = burnup
		use_displaced_mesh = 0
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
    block = 'pellets'
    specific_heat = 320
    thermal_conductivity = 3
  [../]
    [./hcm_sheath]
    type = HeatConductionMaterial_burnup_temp
    block = 'sheath'
    specific_heat = 330
    thermal_conductivity = 16
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

	[./m_density]
    type = Density
		block = 'pellets sheath'
		disp_r = disp_x
		disp_z = disp_y
		density = 10e3
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
		dt = 3
		cutback_factor = 0.5
		growth_factor = 2
		optimal_iterations = 10
		iteration_window = 2 #0.2*optimal_iterations
		linear_iteration_ratio = 100
	[../]

#	[./Quadrature]
#		order = THIRD
#	[../]
[]


[Outputs]
  file_base = rmctests/2D/2D_simple_burnup_out
  output_initial = true
  exodus = true
  print_perf_log = true
[]


#[Debug]
#	show_var_residual_norms = true
#[]
