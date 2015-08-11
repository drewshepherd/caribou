/*SheathThermalMaterial source file
	*
	*Calculates the thermal conductivity, k, for the sheath
	*updates the specific heat at each quadrature point
	*
	*written by Kyle Gamble and Drew Shepherd
*/

#include "SheathThermalMaterial.h"

template<>
InputParameters validParams<SheathThermalMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<bool>("model_thermal_conductivity", "Set true to calculate thermal conductivity");
	params.addRequiredParam<bool>("model_specific_heat", "Set true to calculate specific heat capacity");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
  params.addRequiredCoupledVar("k_sheath", "Coupled thermal conductivity of the sheath");
  params.addRequiredCoupledVar("k_sheath_dT", "Coupled thermal conductivity_dT of the sheath");
  return params;
}

SheathThermalMaterial::SheathThermalMaterial(const InputParameters & parameters) :
    	Material(parameters),
	_model_thermal_conductivity(getParam<bool>("model_thermal_conductivity")),
	_model_specific_heat(getParam<bool>("model_specific_heat")),
  _temp(coupledValue("temp")),
  _grad_temp(coupledGradient("temp")),
  _k(coupledValue("k_sheath")), //the thermal of conductivity for the sheath, exclusively
  _k_dT(coupledValue("k_sheath_dT")),

  _thermal_conductivity(declareProperty<Real>("thermal_conductivity")), //the sheath's contribution to the thermal conductivity of the whole geometry
  _thermal_conductivity_dT(declareProperty<Real>("thermal_conductivity_dT")),
  _specific_heat(declareProperty<Real>("specific_heat"))
{
}

//Compute temperature dependent specific heat of Zircaloy-4
Real 
SheathThermalMaterial::computeSpecificHeat(const Real temp)
{
	Real cp = 0;
	if (temp < 1115)
	{
		cp = 6.55e6 * (1.1061e-4 * temp + 0.2575);
	}	
	else
	{
		cp = 2.3318e6;
	}

	//Divide by the density to get specific heat in J/kg.
	return cp / 6551.0;
}


void
SheathThermalMaterial::computeProperties()
{

	for(_qp=0; _qp<_qrule->n_points(); ++_qp)
  {
		//Thermal conductivity of Zircaloy-4
		_thermal_conductivity[_qp] = _k[_qp];
		_thermal_conductivity_dT[_qp] = _k_dT[_qp];

		//Specific heat of Zircaloy-4 (J/(kg K))
		_specific_heat[_qp] = computeSpecificHeat(_temp[_qp]);
	}
}

	
