###################################################
##2D with heat genration that depends on time 
##taken from heat_source_bar.i
###################################################

[GlobalParams]
	order = SECOND
[]

[Problem]
	type = FEProblem
	coord_type = RZ
[]

[Mesh]
	file = rmctests/2D/2D_simple.e
	displacements = 'disp_x disp_y'
	patch_size = 1000
[]

[Functions]
	[./temperature]
		type = ConstantFunction
		value = 400
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

[Kernels]
	[./heat]
		type = HeatConduction
		variable = temp
	[../]
	[./heat_dt]
		type = HeatConductionTimeDerivative
		variable = temp
#		use_displaced_mesh = 0
	[../]
	[./heatsource]
		type = HeatSource
		block = pellets
		variable = temp
		value = 3e8
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

	[./pelletsz]
    type = DirichletBC
    variable = disp_y
    boundary = bot_land_point
    value = 0
  [../]
  [./pelletsr]
    type = DirichletBC
    variable = disp_x
    boundary = pellets_in
    value = 0
  [../]
	
  [./sheath_z]
    type = DirichletBC
    variable = disp_y
    boundary = 'sheath_bot'
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

[ThermalContact]
  [./pellet_to_sheath]
    slave = pellets_out
    quadrature = true
    master = sheath_in
    variable = temp
    type = GapHeatTransfer
  [../]
[]

[Materials]
	 [./hcm]
    type = HeatConductionMaterial
    block = 'pellets sheath'
    specific_heat = 1000
    thermal_conductivity = 1
  [../]
  [./sm_materials_pellet]
    type = Elastic
		formulation = AxisymmetricRZ
		block = 'pellets'
		youngs_modulus = 100e9
		poissons_ratio = 0.3
		disp_r = disp_x
		disp_z = disp_y
		temp = temp
		thermal_expansion = 1e-5
		stress_free_temperature = 400
  [../]
  [./sm_materials_sheath]
    type = Elastic
		formulation = AxisymmetricRZ
		block = 'sheath'
		youngs_modulus = 10e9
		poissons_ratio = 0.3
		disp_r = disp_x
		disp_z = disp_y
		temp = temp
		thermal_expansion = 0
		stress_free_temperature = 400
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
		variable = temp
		max_increment = 50
		type = MaxIncrement
	[../]
	[./MaxIncrementDispX]
		variable = disp_x
		max_increment = 1e-5
		type = MaxIncrement
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

# controls for linear iterations
  l_max_its = 100
  l_tol = 1e-2

# controls for nonlinear iterations
  nl_max_its = 15
  nl_abs_tol = 1e-10

# time control
#  start_time = 0.0
#  dt = 10.0
#  end_time = 500.0
# num_steps = 5000

	[./TimeStepper]
		type = IterationAdaptiveDT
		
		dt = 0.1
		cutback_factor = 0.5
		growth_factor = 2
		optimal_iterations = 8
		iteration_window = 2 #0.2*optimal_iterations
	[../]

[]

[Outputs]
  file_base = rmctests/2D/2D_simple_heatsourcedt_v1_out
  output_initial = true
  exodus = true
  print_perf_log = true
[]
