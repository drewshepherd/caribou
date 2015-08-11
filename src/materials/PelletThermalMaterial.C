/*PelletThermalMaterial source file
	*
	*Calculates the thermal conductivity, k, the specific heat and porosity at each quadrature point
	*there is no deviation from stoichiometry, k3x = 1
	*
	*written by Kyle Gamble and Drew Shepherd
*/

#include "PelletThermalMaterial.h"

template<>
InputParameters validParams<PelletThermalMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<bool>("model_thermal_conductivity", "Set true to calculate thermal conductivity");
	params.addRequiredParam<bool>("model_specific_heat", "Set true to calculate specific heat capacity");
	params.addRequiredParam<bool>("model_porosity", "Set true to update porosity");
	params.addRequiredParam<bool>("model_alpha", "Set true to model thermal expansion");
	params.addRequiredParam<bool>("model_SFP", "Set true to model SFP");
	params.addRequiredParam<bool>("model_GFP", "Set true to model GFP");
	params.addRequiredParam<bool>("display_values", "display values");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	params.addRequiredCoupledVar("burnup", "Coupled burnup");
	params.addRequiredCoupledVar("burnup_dt", "Coupled burnup_dt");
	params.addRequiredCoupledVar("densification_fraction", "Coupled densification factor");
  params.addRequiredCoupledVar("k_pellets", "Coupled thermal conductivity of the pellets");
  params.addRequiredCoupledVar("k_pellets_dT", "Coupled thermal conductivity_dT of the pellets");
	params.addRequiredParam<Real>("initial_porosity", "The initial porosity");
  return params;
}

PelletThermalMaterial::PelletThermalMaterial(const InputParameters & parameters) 
:Material(parameters),

	_model_thermal_conductivity(getParam<bool>("model_thermal_conductivity")),
	_model_specific_heat(getParam<bool>("model_specific_heat")),
	_model_porosity(getParam<bool>("model_porosity")),
	_model_alpha(getParam<bool>("model_alpha")),
	_model_SFP(getParam<bool>("model_SFP")),
	_model_GFP(getParam<bool>("model_GFP")),
	_display_values(getParam<bool>("display_values")),
  _temp(coupledValue("temp")),
	_burnup(coupledValue("burnup")),
	_burnup_dt(coupledValue("burnup_dt")),
	_densificationF(coupledValue("densification_fraction")),
  _k(coupledValue("k_pellets")), //the pellets thermal conductivity, exclusively
  _k_dT(coupledValue("k_pellets_dT")),

	_initial_porosity(getParam<Real>("initial_porosity")),

	_density(getMaterialProperty<Real>("density")),
	_enrichment(getMaterialProperty<Real>("enrich")),
	_thermal_conductivity(declareProperty<Real>("thermal_conductivity")), //the pellets's contribution to the thermal conductivity of the whole geometry
  _thermal_conductivity_dT(declareProperty<Real>("thermal_conductivity_dT")),
  _specific_heat(declareProperty<Real>("specific_heat")),
  _alpha1(declareProperty<Real>("alpha")),
  _alpha1_old(declarePropertyOld<Real>("alpha")),
  _SFP(declareProperty<Real>("solid_products")),
  _GFP(declareProperty<Real>("gaseous_products")),
  _porosity(declareProperty<Real>("porosity"))
{
}

//defines the initial values
void
PelletThermalMaterial::initStatefulProperties(unsigned n_points)
{
  for (unsigned qp(0); qp < n_points; ++qp)
  {
    _porosity[qp] = _initial_porosity;
  }
}

//compute the theoretical density
Real 
PelletThermalMaterial::computeTheoDensity(const Real alpha, const Real temp)
{
	const Real rho_273 = 10963.0;
	
	const Real rho_th = rho_273 * std::pow(alpha * (temp - 273) + 1, -3);

	return rho_th;
}

//Compute Temperature Dependent Specific Heat
Real
PelletThermalMaterial::computeSpecificHeat(const Real temp)
{
//many of these values can be found in Table 7 of Prudil
	const Real R = 1.987; 				//cal/(mol K)
	const Real theta = 535.285;		// K, 548.68 J/(kg K)
	const Real Ed = 37694.6;			//cal/mol
	const Real k1 = 19.145;				//cal/(mol K), 302.27 J/(kg K)
	const Real k2 = 7.8473e-4;		//cal/(mol K^2), 8.463e-3 J/(kg K^2)
	const Real k3 = 5.6437e6;			//cal/mol, 8.741e7 J/kg

	const Real temp2 = temp * temp;

	const Real cp = (k1*std::pow(theta,2)*std::exp(theta/temp))/(temp2*std::pow(std::exp(theta/temp)-1,2) + 2*k2*temp + ((k3*Ed)/(R*temp2))*std::exp(-Ed/(R*temp)));

	//Convert to SI units (J/kg)
	return cp * 15.496;
}


//calculates alpha
Real
PelletThermalMaterial::computeAlpha(const Real temp)
{
	const Real K1 = 1.0e-5;				//K^-1
	const Real K2 = 3.0e-3;
	const Real K3 = 4.0e-2;
	const Real k = 1.3806e-23;		//Boltzmann Constant (J/K)
	const Real Ed = 6.9e-20;			//J
	
	const Real temp2 = temp * temp;
	const Real alpha = (K1 + ((K3 * Ed) / (k * temp2)) * std::exp(-Ed / (k * temp))) / (1 + K1 * temp - K2 + K3 * std::exp(-Ed / (k * temp)));
	return alpha;
}

//calculates SFP
Real
PelletThermalMaterial::computeSFP(const Real burnup)
{
	return 0.0032 * burnup / 225; //eq: 5.56 from Prudil
}

//calculates GFP
Real
PelletThermalMaterial::computeGFP(const Real density, const Real enrichment, const Real burnup, const Real burnup_dt, const Real temp)
{
	const Real MU235 = 235.0439; //mass of U-235
	const Real MU238 = 238.0508; //mass of U-238
	const Real MUO2 = 270.03; //mass of UO2

	const Real density_u = density * (MU235 * enrichment + MU238 * (1 - enrichment)) / MUO2; //density of uranium
	const Real E_f = 200 * (1e6 * 1.602e-19); //MeV, energy per fission in [J]
	const Real Bu_f = burnup * 3.6e3 * 1e6 * density_u / E_f; //number of fissions per cubic metre
	const Real Bu_f_dt = burnup_dt * 3.6e3 * 1e6 * density_u / E_f; //number of fissions per cubic metre per time

	const Real value = std::pow((2800 - temp), 11.73) * exp(-0.0162 * (2800 - temp)- 8e-27 * Bu_f) * Bu_f_dt * _dt; //eq: 5.57 from prudil

	return 8.8e-56 * value;
}

//Compute Properties
void
PelletThermalMaterial::computeProperties()
{
	for(_qp=0; _qp<_qrule->n_points(); ++_qp)
  {
		const Real temp_ = _temp[_qp];
		const Real burnup_ = _burnup[_qp];
		const Real burnup_dt_ = _burnup_dt[_qp];
		const Real densificationF_ = _densificationF[_qp];
		const Real enrichment_ = _enrichment[_qp];
		const Real density_ = _density[_qp];

//modelling thermal expansion coefficient
		if (_model_alpha)
		{
			_alpha1[_qp] = computeAlpha(temp_);
		}
		else if (!_model_alpha)
		{
			_alpha1[_qp] = 1;//giving default value of 1
		}

//modelling solid fission products swelling strain
		if (_model_SFP)
		{
			_SFP[_qp] = computeSFP(burnup_);
		}
		else if (!_model_SFP)
		{
			_SFP[_qp] = 1;//giving default value of 1
		}

//modelling gaseous fission products swelling strain
		if (_model_GFP)
		{
			_GFP[_qp] = computeGFP(density_, enrichment_, burnup_, burnup_dt_, temp_);
		}
		else if (!_model_GFP)
		{
			_GFP[_qp] = 1;//giving default value of 1
		}

//Density
		Real rho_theo = 0.0;

		if (_t_step == 0 || !_model_porosity)
		{
			_porosity[_qp] = _initial_porosity;
		}
		else if (_t_step != 0 && _model_porosity && _model_GFP)
		{
			rho_theo = computeTheoDensity(_alpha1[_qp], temp_);
			_porosity[_qp] = _initial_porosity * (1 - densificationF_) + _GFP[_qp];
		}
		else if (_t_step != 0 && _model_porosity && !_model_GFP)
		{
			rho_theo = computeTheoDensity(_alpha1[_qp], temp_);
			_porosity[_qp] = _initial_porosity * (1 - densificationF_);
		}


//modelling thermal conductivity
		if (_model_thermal_conductivity)
		{
			_thermal_conductivity[_qp] = _k[_qp];
			_thermal_conductivity_dT[_qp] = _k_dT[_qp];
		}
		else if (!_model_thermal_conductivity)
		{
			_thermal_conductivity[_qp] = 1; //giving default value of 1
			_thermal_conductivity_dT[_qp] = 1; //giving default value of 1
		}

//modelling specific heat
		if (_model_specific_heat)
		{
			//Specific heat of UO2 (J/(kg K))
			_specific_heat[_qp] = computeSpecificHeat(temp_);
		}
		else if (!_model_specific_heat)
		{
			_specific_heat[_qp] = 1; //giving default value of 1
		}

		if (_display_values && _burnup[_qp] > 1e-4)
		{
		 	std::stringstream msg;
    	msg << "\tspecific heat: " << _specific_heat[_qp] << "\n"
					<< "\tthermal cond: " << _thermal_conductivity[_qp] << "\n"
					<< "\talpha: " << _alpha1[_qp] << "\n"
					<< "\talpha_old: " << _alpha1_old[_qp] << "\n"
					<< "\talpha_increment: " << _alpha1[_qp] - _alpha1_old[_qp] << "\n"
					<< "\tsfp: " << _SFP[_qp] << "\n"
					<< "\tgfp: " << _GFP[_qp] << "\n"
					<< "\tburnup: " << burnup_ << "\n"
					<< "\tdensification: " << densificationF_ << "\n"
					<< "\tburnup_dt: " << burnup_dt_ << "\n"
					<< "\tporosity: " << _porosity[_qp] << "\n";
    	mooseWarning( msg.str() );
		}
	}
}
