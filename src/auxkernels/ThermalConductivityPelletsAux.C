/*ThermalConductivityPellets AuxKernel source file
	*
	*A method of calculating thermal conductivity
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "ThermalConductivityPelletsAux.h"


template<>
InputParameters validParams<ThermalConductivityPelletsAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_k", "Set true to calculate Thermal conductivity of the sheath");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	params.addRequiredCoupledVar("burnup", "Coupled burnup");
	return params;
}

ThermalConductivityPelletsAux::ThermalConductivityPelletsAux(const InputParameters & parameters)
  :AuxKernel(parameters),

	_model_k(getParam<bool>("model_k")),
  _temp(coupledValue("temp")),
	_burnup(coupledValue("burnup")),
	
	_porosity(getMaterialProperty<Real>("porosity"))
{
}

Real
ThermalConductivityPelletsAux::computeValue()
{
	if (_model_k)
	{
		const Real temp_ = _temp[_qp];
		const Real burnup_ = _burnup[_qp];
		const Real porosity_ = _porosity[_qp];

		//Thermal conductivity of UO2
		Real cond0;
		Real kappa1d;
		Real kappa1p;
		Real kappa2p;
		Real kappa4r;

		computeChassieUnirradThCond(temp_, cond0);
		computeDissFissionProdThCond(temp_, burnup_, kappa1d);
		computePrecFissionProdThCond(temp_, burnup_, kappa1p);
		computePoresThCond(temp_, porosity_, kappa2p);
		computeRadDamageThCond(temp_, kappa4r);

		const Real k = cond0 * kappa1d * kappa1p * kappa2p * kappa4r;

		return k;
	}
	else
	{
		return 1; //default value
	}

}

//Compute Unirradiated Thermal Conductivity
void
ThermalConductivityPelletsAux::computeChassieUnirradThCond(const Real temp, Real & cond0)
{
	const Real A = 0.030771;			//(m K)/W
	const Real B = 2.25e-4;				//m/W
	const Real C = 9.28e9;				//(W K)/m
	const Real D = 18295.09;			//K

	const Real temp2 = temp * temp;
	
	cond0 = (1 /(A + B * temp)) + (C / temp2) * std::exp(-D / temp);
}

//Compute contribution to thermal conductivity due to dissolved fission products
void
ThermalConductivityPelletsAux::computeDissFissionProdThCond(const Real temp, const Real burnup, Real & kappa1d)
{
	Real n = 0.0;
	Real m = 0.0;
	Real beta = burnup / 225;	//convert to atom percent

	if (beta <=0.1)
	{	
		kappa1d = 1.0;
	}
	else
	{
		beta = burnup / 225.;	//convert to atom percent
		m = 1.09 / std::pow(beta, 3.265);
		n = 0.0643 / std::sqrt(beta);
	
		kappa1d = (m + n * std::sqrt(temp)) * std::atan(1 / (m + n * std::sqrt(temp)));
	}
}

//Compute contribution to thermal conductivity due to solid fission products
void
ThermalConductivityPelletsAux::computePrecFissionProdThCond(const Real temp, const Real burnup, Real & kappa1p)
{
	Real beta = burnup / 225;  //convert to atom percent
	
	kappa1p = 1 + ((0.019 * beta) / (3 - 0.019 * beta)) * (1 / (1 + std::exp(-(temp - 1200.) / 100.)));
}

//Computer contribution to thermal conductivity due to pores and fission gas bubbles
void
ThermalConductivityPelletsAux::computePoresThCond(const Real temp, const Real porosity, Real & kappa2p)
{
	kappa2p = 1.0 - (2.05 - 5.0e-4 * temp) * porosity;
}

//Computer contribution to thermal conductivity due to radiation damage
void
ThermalConductivityPelletsAux::computeRadDamageThCond(const Real temp, Real & kappa4r)
{
	kappa4r = 1 - 0.2 / (1 + std::exp((temp - 900.) / 80.));
}

