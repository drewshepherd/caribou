###################################################
##3D with Q_fission heat generation, porosity and burnup-dependent heat conduction, thermal expansion OFF
###################################################
##***********contact causes a segmentation fault

[GlobalParams]
	family = LAGRANGE
	order = FIRST
	#disp_x = disp_x
	#disp_y = disp_y
[]

##Trelis geometry
[Mesh]
	file = 3D_Geometries/3D_endcapcoarse_v1.e #endcap_v1.e
	displacements = 'disp_x disp_y disp_z'
	patch_size = 1000
[]

##Functions to be used in the Kernels, AuxKernels and Materials
[Functions]
	[./linear_power_function]
		type = PiecewiseLinear
		x = '0 10'
		y = '0 43000'
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
[]

##Variables to be solved for in the AuxKernels
[AuxVariables]
	[./avg_burnup] ##solves for the average burnup - only used in FissionHeatMaterial
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./burnup] ##solves for the local burnup - used everywhere except FissionHeatMaterial
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./burnup_dt] ##solves for the local burnup_dt
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./vstrain] ##solves for the volumetric strain
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./k_pellets] 
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./k_pellets_dT] 
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
	[./k_sheath] 
		order = FIRST 
		family = MONOMIAL
		block = 'sheath endcap'
	[../]
	[./k_sheath_dT] 
		order = FIRST 
		family = MONOMIAL
		block = 'sheath endcap'
	[../]
	[./test]
		order = FIRST 
		family = MONOMIAL
		block = pellets
	[../]
[]
##Variables to be solved for in the differential equations (ODE/PDE) in the Kernels
[Variables]
	[./temp] ##temperature variable
		initial_condition = 550
		block = 'pellets sheath endcap'
	[../]
	[./disp_x] ##x-displacement variable
  [../]
  [./disp_y] ##y-displacement variable
  [../]
  [./disp_z] ##z-displacement variable
  [../]
	[./densification_fraction] ##densification fraction variable - a percentage less than 60%
		block = pellets
		initial_condition = 0
	[../]
[]

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
		initial_porosity = 0.0286
		execute_on = timestep_begin
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
		model_density = false
		model_q_fission = false
		model_q_fission_old = false
		model_ratio = false
		model_pellet_radius = false
		model_porosity = false
		model_enrichment = false
		model_SFP = false
		model_GFP = true
	[../]
[]

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
		block = pellets
	[../]
	[./Densification_Fraction_RHS] ##solves for F_porosity - couples to Denisification_Fraction_LHS
		type = DensificationFractionKernel
		variable = densification_fraction
		temp =temp
		burnup_dt = burnup_dt
		model_densification_fraction = true
	[../]
	[./Densification_Fraction_LHS] ##solves for F_porosity - couples to Denisification_Fraction_RHS
		type = TimeDerivative
		variable = densification_fraction
	[../]
[]

##Allows for solid mechanics to be used
[SolidMechanics]
  [./solid]
	  temp = temp
  	disp_x = disp_x
		disp_y = disp_y
  	disp_z = disp_z		 
  [../]
[]


[BCs]
	[./Sheath_out] ##Specifies the coolant temperature and the rate of heat flux out of the bundle
		type = ConvectiveFluxFunction
    variable = temp
    boundary = 'sheath_outer endcap_out endcap_top'
		T_infinity = 550 
		coefficient = 5e4
	[../]
	[./Z_up] ##Constrains the top of the pellets and sheath to not move vertically
		type = DirichletBC
		variable = disp_z
		boundary = 'pellet1_bot sheath_bot'
		value = 0
	[../]
	[./X_in] ##Constrains the x-halfplane of the pellet to not move horizontally
		type = DirichletBC
		variable = disp_x
		boundary = 'pellets_in sheath_in endcap_in'
		value = 0
	[../]
	[./Y] ##Constrains the y movement
		type = DirichletBC
		variable = disp_y
		boundary = 'pellet1_point'
		value = 0
	[../]
#	[./Pressure]
#		[./coolant_pressure]
#			disp_x = disp_x
#			disp_y = disp_y
#  		disp_z = disp_z	
#			boundary = sheath_out
#			factor = 1e7
#			function = pressure_ramp
#  	[../]
#	[../]
[]

##Specifies the properties when the pellet-sheath and pellet-endcap contact each other - specifies the penalty factor coefficient
#[Contact]
#	[./pellet_to_sheath_contact]
#		slave = pellets_out
#   master = sheath_IN
#		formulation = penalty
#		model = frictionless
#		penalty = 1e8
# 	disp_x = disp_x
#		disp_y = disp_y
#  	disp_z = disp_z		
#	[../]
#	[./pellet_to_endcap_contact]
#		slave = pellet1_bot
#    master = endcap_top
#		formulation = penalty
#		model = frictionless
#		penalty = 1e10
# 	disp_x = disp_x
#		disp_y = disp_y
#  	disp_z = disp_z		
#  [../]
#[]

##Allows heat to transfer over the gap between the pellet-sheath, pellet-endcap, and sheath-endcap
[ThermalContact]
  [./pellet_to_sheath]
		type = GapHeatTransfer
    variable = temp
    slave = pellets_out
    master = sheath_inner
		quadrature = true
		gap_conductivity= 0.15
  [../]
	[./pellet_to_endcap]
		type = GapHeatTransfer
    variable = temp
    slave = pellet2_top
    master = endcap_bot
		quadrature = true
  [../]
  [./pellet_to_pellet]
		type = GapHeatTransfer
		variable = temp
		slave = pellet2_bot
    master = pellet1_top
		quadrature = true
  [../]
[]


[Materials]
	[./HCM_Pellet] ##Determines the thermal properties of the pellet - thermal conductivity, specific heat, porosity. If these models are set to false, then default values of 1 are given
    type = PelletThermalMaterial
		temp = temp
		burnup = burnup
		burnup_dt = burnup_dt
		densification_fraction = densification_fraction
		k_pellets = k_pellets
		k_pellets_dT = k_pellets_dT
		block = pellets
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
  [./HCM_Endcap] ##Determines the thermal properties of the endcap - thermal conductivity, specific heat, thermal expansion
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
		block = pellets
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
	[./Mechanical_Endcap] ##Determines the solid mechanical properties of the endcap - stresses, strains
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
		block = pellets
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
	[./Endcap_density] ##Specifies the density of the endcap
    type = Density
		block = endcap
	 	disp_x = disp_x
		disp_y = disp_y
  	disp_z = disp_z
		density = 6.5e3
  [../]
	[./FissionHeat] ##Determines the amount of heat produced by nuclear fission
		type = FissionHeatMaterial
		burnup_avg = avg_burnup
		block = pellets
		linear_power = linear_power_function_time #burnup
		initial_fuel_density = 1.065e4
		initial_qfission = 0
		enrichment = 1
		pellet_radius = 6e-3
		ratio = 0.925
		is_3D = true
		model_Qfission = true
[]


##Specifies the maximum increment each of these variables may take in each time step
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
		max_increment = 3e-6
		type = MaxIncrement
	[../]
	[./MaxIncrementDensificationF]
		variable = densification_fraction
		max_increment = 1e-2
		type = MaxIncrement
	[../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
   full = true
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
  nl_abs_tol = 5e-10 #10
	nl_rel_tol = 1e-6 #5

#time control
  start_time = 0.0
#  dt = 10.0
  end_time = 1e8
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
  file_base = 3D_Paraview/3D_endcap_v1_out
  output_initial = true
  exodus = true
  print_perf_log = true
[]


#[Debug]
#	show_var_residual_norms = true
#[]
