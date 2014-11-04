/*StrainMaterial source file
	*
	*Calculates stress, strain, young's modulus and poisson's ratio for the pellets
	*Calculates the thermal expansion coefficient to allow for thermal deformation (strain) - this is done by updating the "modifystrainincrement" built into MOOSE
	*
	*written by Kyle Gamble and Drew Shepherd
*/

#include "StrainMaterial.h"
#include "SymmIsotropicElasticityTensor.h"

template<>
InputParameters validParams<StrainMaterial>()
{
  InputParameters params = validParams<SolidModel>();

  params.addRequiredParam<bool>("model_thermal_expansion", "Set true to turn on thermal expansion model");
	params.addRequiredParam<bool>("model_youngs_modulus", "Set true to calculate elastic moduli internally");
  params.addRequiredParam<bool>("display_values", "display values");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	params.addRequiredCoupledVar("burnup", "Coupled burnup");
	params.addRequiredCoupledVar("burnup_dt", "Coupled burnup_dt");
	params.addRequiredCoupledVar("vstrain", "Coupled volumetric strain");
  return params;
}

StrainMaterial::StrainMaterial( const std::string & name, InputParameters parameters ) :
  SolidModel( name, parameters ),

  _model_thermal_expansion(getParam<bool>("model_thermal_expansion")),
	_model_youngs_modulus(getParam<bool>("model_youngs_modulus")),
	_display_values(getParam<bool>("display_values")),

  _temp(coupledValue("temp")),
  _temp_old(coupledValueOld("temp")),
	_burnup(coupledValue("burnup")),
	_burnup_old(coupledValueOld("burnup")),
	_burnup_dt(coupledValue("burnup_dt")),
	_burnup_dt_old(coupledValueOld("burnup_dt")),
	_vstrain(coupledValue("vstrain")),
	_vstrain_old(coupledValueOld("vstrain")),

  _thermal_strain_increment(0.0),
	_alpha1(getMaterialProperty<Real>("alpha")),
	_alpha1_old(getMaterialPropertyOld<Real>("alpha")),
	_density(getMaterialProperty<Real>("density")),
	_density_old(getMaterialPropertyOld<Real>("density")),
	_enrichment(getMaterialProperty<Real>("enrich"))
{
}

void
StrainMaterial::computeStress()
{
	if(_t_step == 0) return;
	SymmTensor stress_new( (*elasticityTensor()) * _strain_increment );
	_stress[_qp] = stress_new;
	_stress[_qp] += _stress_old;
}

Real
StrainMaterial::computeYoungsModulus(const Real temp)
{	
	const Real D = _density[_qp] / 10963.0;

	const Real YM = 2.334e11 * (1 - 2.752 * (1 - D)) * (1 - 1.0915e-4 * temp);
	
	return YM;
}

void
StrainMaterial::modifyStrainIncrement()
{
	if( _t_step != 0 && _has_temp &&  _model_thermal_expansion)
	{
		const Real temp(_temp[_qp]);
		const Real temp_old(_temp_old[_qp]);
  	const Real alpha(_alpha1[_qp]);
  	const Real alpha_old(_alpha1_old[_qp]);

  	const Real alpha_avg = (alpha + alpha_old) / 2.0; //calculated for the most recent time interval
		_thermal_strain_increment  = alpha_avg * (temp - temp_old);

 		const Real dens_dT(0.0);
		const Real GFP_dT = computeGFP_dT(_density[_qp], _enrichment[_qp], _burnup[_qp], _burnup_dt[_qp], temp); //calculated for the most recent time interval

  	const Real dVStrain_dT(dens_dT + GFP_dT);

		const Real vstrain_increment = _vstrain[_qp] - _vstrain_old[_qp];
		_strain_increment.addDiag(-vstrain_increment / 3 - _thermal_strain_increment);

		_d_strain_dT.zero();
   	_d_strain_dT.addDiag( -dVStrain_dT / 3 - alpha_avg);

		if (_display_values && _burnup[_qp] > 1e-4)
		{			
		 	std::stringstream msg;
     	msg << "\tstrain incremement: " << _strain_increment << "\n"
					<< "\tthermal_strain: " << _thermal_strain_increment << "\n"
					<< "\tvstrain increment: " << vstrain_increment << "\n"
					<< "\talphaincrement: " << alpha - alpha_old << "\n"
					<< "\ttempincrement: " << temp - temp_old << "\n"
					<< "\tburnup increment: " << _burnup[_qp] - _burnup_old[_qp] << "\n"
					<< "\tburnup _dtincrement: " << _burnup_dt[_qp] - _burnup_dt_old[_qp] << "\n"
					<< "\tdensity increment: " << _density[_qp] - _density_old[_qp] << "\n"
					<< "\t_d_strain_dT: " << _d_strain_dT << "\n";
     	mooseWarning( msg.str() );
		}
	}
}

Real
StrainMaterial::computeGFP_dT(const Real density, const Real enrichment, const Real burnup, const Real burnup_dt, const Real temp)
{
	const Real MU235 = 235.0439; //mass of U-235
	const Real MU238 = 238.0508; //mass of U-238
	const Real MUO2 = 270.03; //mass of UO2

	const Real density_u = density * (MU235 * enrichment + MU238 * (1 - enrichment)) / MUO2; //density of uranium
	const Real E_f = 200 * (1e6 * 1.602e-19); //MeV, energy per fission in [J]
	const Real Bu_f = burnup * 3.6e3 * 1e6 * density_u / E_f; //number of fissions per cubic metre
	const Real Bu_f_dt = burnup_dt * 3.6e3 * 1e6 * density_u / E_f; //number of fissions per cubic metre per time

	const Real value = Bu_f_dt * std::exp(-0.0162 * (2800 - temp) * Bu_f) * (0.0162 * Bu_f * std::pow(2800 - temp, 11.73) - 11.73 * std::pow(2800 - temp, 10.73));

	return 8.8e-56 * value;
}

bool
StrainMaterial::updateElasticityTensor(SymmElasticityTensor & tensor)
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
    const Real PR = 0.316;

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
