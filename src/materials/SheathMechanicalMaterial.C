/*PelletMechanicalMaterial source file
	*
	*Calculates stress, strain, young's modulus and poisson's ratio for the sheath
	*Calculates the thermal expansion coefficient to allow for thermal deformation (strain)
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "SheathMechanicalMaterial.h"
#include "SymmIsotropicElasticityTensor.h"

template<>
InputParameters validParams<SheathMechanicalMaterial>()
{
   InputParameters params = validParams<SolidModel>();

   params.addRequiredCoupledVar("temp", "coupled temperature");
	//Sub-Newton Iteration control parameters
   params.addRequiredParam<Real>("relative_tolerance", "Relative convergence tolerance for sub-newtion iteration");
   params.addRequiredParam<Real>("absolute_tolerance", "Absolute convergence tolerance for sub-newtion iteration");
   params.addRequiredParam<unsigned int>("max_its", "Maximum number of sub-newton iterations");
   params.addRequiredParam<bool>("output_iteration_info", "Set true to output sub-newton iteration information");

  //Model options
  params.addRequiredParam<bool>("model_thermal_expansion", "Set true to turn on thermal expansion model");
	params.addRequiredParam<bool>("model_youngs_modulus", "Set true to calculate elastic moduli internally");
	params.addRequiredParam<bool>("model_diffusional_creep", "Set true to turn on diffusional creep model");
  params.addRequiredParam<bool>("display_values", "display values");
   return params;
}

SheathMechanicalMaterial::SheathMechanicalMaterial( const std::string & name,
                                                        InputParameters parameters ) :
  SolidModel( name, parameters ),
	_temp(coupledValue("temp")),
	_relative_tolerance(parameters.get<Real>("relative_tolerance")),
	_absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
	_max_its(parameters.get<unsigned int>("max_its")),
	_output_iteration_info(getParam<bool>("output_iteration_info")),
  _model_thermal_expansion(getParam<bool>("model_thermal_expansion")),
	_model_youngs_modulus(getParam<bool>("model_youngs_modulus")),
	_model_diffusional_creep(getParam<bool>("model_diffusional_creep")),
	_display_values(getParam<bool>("display_values")),

	_creep_strain(declareProperty<SymmTensor>("creep_strain")),
	_creep_strain_old(declarePropertyOld<SymmTensor>("creep_strain")),
	_primary_creep_strain(declareProperty<Real>("primary_creep_strain")),
	_primary_creep_strain_old(declarePropertyOld<Real>("primary_creep_strain"))

{
}

void
SheathMechanicalMaterial::computeStress()
{
	if(_t_step == 0) return;

	// compute trial stress
	SymmTensor stress_new( (*elasticityTensor()) * _strain_increment );
	stress_new += _stress_old;

	// compute deviatoric trial stress
	SymmTensor dev_trial_stress(stress_new);
	dev_trial_stress.addDiag( -stress_new.trace()/3.0);

	// compute effective trial stress
	Real dts_squared = dev_trial_stress.doubleContraction(dev_trial_stress);
	Real effective_trial_stress = std::sqrt(1.5*dts_squared);

	// Use Newton sub-iteration to determine effective creep strain increment

  unsigned int it = 0;
  Real creep_residual = 10.;

  Real del_p(0);
  Real norm_residual = 10.;
  Real first_norm_residual = 10.;
	Real creeprate(0.0);
  Real dphi_ddelp(0);
	Real edotgb(0.0);

	Real T = _temperature[_qp];

  while(it < _max_its
        && norm_residual > _absolute_tolerance
        && (norm_residual/first_norm_residual) > _relative_tolerance )

  {
    const Real stress_delta = effective_trial_stress - 3.*_shear_modulus*del_p;
    creeprate = 0.0;
    dphi_ddelp = 0.0;
	
		//Diffusional Creep
		if(_model_diffusional_creep)
		{
			creeprate += EDOTGB(T, stress_delta/1.0e6);
			dphi_ddelp += -3.0 * _shear_modulus/1.0e6 * d_EDOTGB_d_sigma(T, stress_delta / 1.0e6);
		}

		creep_residual = creeprate -  del_p / _dt;
    norm_residual  = std::abs(creep_residual);

		if(it==0) first_norm_residual = norm_residual;

    del_p = del_p + (creep_residual / (1 / _dt - dphi_ddelp));
	
		//iteration output
		if (_output_iteration_info)
    {
      std::cout
        <<" it=" <<it
        <<" dt=" <<_dt
        <<" temperature=" << _temperature[_qp]
        <<" trial stress=" <<effective_trial_stress
        <<" creeprate=" <<creeprate
        <<" dphi=" <<dphi_ddelp
        <<" creep_residual=" <<norm_residual
        <<" del_p=" <<del_p
        <<" relative tolerance=" << _relative_tolerance
        <<" absolute tolerance=" << _absolute_tolerance
				<<" strain increment = " << _strain_increment
        <<"\n";
    }

    it++;
  }
	if(it == _max_its && (norm_residual/first_norm_residual) > _relative_tolerance && norm_residual > _absolute_tolerance)
  {
    std::cerr
      << " it = " << it << std::endl
      << " dt = " << _dt << std::endl
      << " temperature = " << _temperature[_qp] << std::endl
      << " trial stress = " << effective_trial_stress << std::endl
      << " creeprate = " << creeprate << std::endl
      << " dphi = " << dphi_ddelp << std::endl
      << " creep_residual = " << norm_residual << std::endl
      << " del_p = " << del_p << std::endl
      << " relative tolerance=" << _relative_tolerance
      << " absolute tolerance=" << _absolute_tolerance
      << std::endl;
    mooseError("Max sub-newton iteration hit during creep solve!");
  }

// compute creep and elastic strain increments (avoid potential divide by zero - how should this be done)?
  if (effective_trial_stress < 0.01)
  {
    effective_trial_stress = 0.01;
  }
  SymmTensor creep_strain_increment(dev_trial_stress);

  creep_strain_increment *= (1.5 * del_p / effective_trial_stress);

  _strain_increment -= creep_strain_increment;
	// update stress and creep strain
// compute stress increment

  _stress[_qp] =  *elasticityTensor() * _strain_increment;
  _stress[_qp] += _stress_old;

  _creep_strain[_qp] = creep_strain_increment;
  _creep_strain[_qp] += _creep_strain_old[_qp];
}	

double
SheathMechanicalMaterial::EDOTGB(double T, double sigma_a)
{
	double G = (1000 * (36.3 - 0.0223 * (T - 273)));
	double F = 6.34e6 / (std::pow(G, 2));
	const double R = 8.3144621;
	double m = 0.0;
	double d = 0.0;
	double Q = 0.0;
	if (T <= 1073)  															//Alpha Phase
	{ 
		m = 2.0;
		d = 3;																//must be in micrometers, not meters
		Q = 9431;
	}
	else if ((T > 1073) && (T <= 1273))						//Transition Phase
	{
		m = -5e-4*T + 2.5365;
		d = 0.485*T - 517.405;								//must be in micrometers, not meters
		Q = -16.96*T + 27629.08;
	}
	else 																					//Beta Phase
	{
		m = 1.9;
		d = 100;															//must be in micrometers, not meters
		Q = 6039;
	}
		
	return F * std::pow(sigma_a/d, m)*std::exp(-Q/T);
}

double
SheathMechanicalMaterial::d_EDOTGB_d_sigma(double T, double sigma_a)
{
	double G = (1000 * (36.3 - 0.0223 * (T - 273)));
	double F = 6.34e6/(std::pow(G, 2));
	const double R = 8.3144621;
	double m = 0.0;
	double d = 0.0;
	double Q = 0.0;

	if (T <= 1073)  															//Alpha Phase
	{ 
		m = 2.0;
		d = 3;																//must be in micrometers, not meters
		Q = 9431;
	}
	else if ((T > 1073) && (T <= 1273))						//Transition Phase
	{
		m = -5e-4 * T + 2.5365;
		d = 0.485 * T - 517.405;								//must be in micrometers, not meters
		Q = -16.96 * T + 27629.08;
	}
	else 																					//Beta Phase
	{
		m = 1.9;
		d = 100;															//must be in micrometers, not meters
		Q = 6039;
	}
		
	return F * m * std::pow(sigma_a, m - 1) * std::pow(1 / d, m) * std::exp(-Q / T);
}		

Real
SheathMechanicalMaterial::computeAxialThermEx(const Real temp)
{
	Real Athex = 0;
	if (temp <= 1073) 													//Alpha Phase
	{
		Athex = -2.506e-5 + 4.441e-6 * (temp - 273); //eq: 6.61 from Prudil
	}
	else if ((temp > 1073) && (temp <= 1273))  //Transition Phase
	{
		Athex = 0.0120387 - 1.06387e-5 * (temp - 273); //eq: 6.65 from Prudil
	}
	else 																			//Beta Phase
	{
		Athex = -8.3e-3 + 9.7e-6 * (temp - 273);	//eq: 6.63 from Prudil
	}
	return Athex;
}

Real
SheathMechanicalMaterial::computeRadialThermEx(const Real temp)
{
	Real Dthex = 0;
	if (temp <= 1073)
	{
		Dthex = -2.373e-4 + 6.721e-6 * (temp - 273); //eq: 6.62 from Prudil
	}
	else if ((temp > 1073) && (temp <= 1273))
	{
		Dthex = 0.0140975 - 1.11975e-5 * (temp - 273); //eq: 6.66 from Prudil
	}
	else
	{
		Dthex = -6.8e-3 + 9.7e-6 * (temp - 273); //eq: 6.64 from Prudil
	}

	return Dthex;
}

Real
SheathMechanicalMaterial::computeYoungsModulus(const Real temp)
{
	Real YM = 0;
	if (temp < 1135)
	{
		YM = 1.148e11 - 5.99e7 * temp;
	}
	else if ((temp >= 1135) && (temp < 2120))
	{
		YM = 1.005e11 - 4.725e7 * temp;
	}
	else
	{
		YM = 3.30e8;
	}
	
	return YM; //eq: 6.60 from Prudil
}

void
SheathMechanicalMaterial::modifyStrainIncrement()
{
  if (_model_thermal_expansion && _t_step != 0)
  {
			const Real temp(_temperature[_qp]);
      Real temp0(_temperature_old[_qp]);

			const Real Athex = computeAxialThermEx(temp);
			const Real Athex0 = computeAxialThermEx(temp0);
			const Real Dthex = computeRadialThermEx(temp);
			const Real Dthex0 = computeRadialThermEx(temp0);

		  SymmTensor thermal_strain_increment;
		  thermal_strain_increment.zero();
			//the axial and radial strains are different (the sheath is non-isotropic)
		  thermal_strain_increment(0,0) = Dthex - Dthex0;
		  thermal_strain_increment(1,1) = Dthex - Dthex0;
		  thermal_strain_increment(2,2) = Athex - Athex0;

		  _strain_increment -= thermal_strain_increment;

			if (_display_values)// && _temp[_qp] > 1e3)
			{
		 		std::stringstream msg;
     		msg << "\tstrain incremement: " << _strain_increment << "\n"
						<< "\ttemp: " << temp << "\n"
						<< "\ttemp old: " << temp0 << "\n"
						<< "\ttemp increment: " << temp - temp0 << "\n";
     		mooseWarning( msg.str() );
			}
	}
	else
  {
    if ( _t_step != 0 )
    {
       const Real tStrain( _alpha * (_temperature[_qp] - _temperature_old[_qp]) ); //where does alpha come from?
       _strain_increment.addDiag( -tStrain );

       _d_strain_dT.zero();
       _d_strain_dT.addDiag( -_alpha );
    }
  }
}

bool
SheathMechanicalMaterial::updateElasticityTensor(SymmElasticityTensor & tensor)
{
  if( _model_youngs_modulus )
  {
		const Real temp = _temperature[_qp];

    SymmIsotropicElasticityTensor * t = dynamic_cast<SymmIsotropicElasticityTensor*>(&tensor);
    if (!t)
    {
      mooseError("Cannot use Youngs modulus or Poissons ratio functions");
    }
    t->unsetConstants();

    const Real YM = computeYoungsModulus(temp);
    const Real PR = 0.3;

    t->constant(false);
    t->setYoungsModulus(YM);
    t->setPoissonsRatio(PR);

    bool changed = true;
    return changed;
  }
  else
  {
    return SolidModel::updateElasticityTensor( tensor );
  }
}


