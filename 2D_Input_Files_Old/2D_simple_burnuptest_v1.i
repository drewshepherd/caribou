###################################################
##2D with heat genration that depends on time 
##taken from heat_source_bar.i
###################################################

[GlobalParams]
	family = LAGRANGE
	order = SECOND
[]


[Problem]
	coord_type = RZ
[]


[Mesh]
	file = rmctests/2D/2D_simple_v1.e
[]


[Functions]
	[./temperature]
		type = ConstantFunction
		value = 400
	[../]
	[./body]
		type = ParsedFunction
		value = 10*x #there is no heat generation at sheath_out, and maximum generation at pellets_in
		#value = 3e8
	[../]
[]


[Variables]
	[./u]
		initial_condition = 0
		block = pellets
	[../]
[]


[Kernels]
	[./BurnupSource] ## [MWh/kgU], fraction of U atoms that have been fissioned
		type = BodyForce
		variable = u
		block = 'pellets'
		function = body
		value = 1
	[../]

	[./burnupTime]
		type = TimeDerivative
		variable = u
		block = 'pellets'
	[../]
[]


[BCs]
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
  file_base = rmctests/2D/2D_simple_burnuptest_v1_out
  output_initial = true
  exodus = true
  print_perf_log = true
[]


#[Debug]
#	show_var_residual_norms = true
#[]
