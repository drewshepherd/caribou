/*ThermalConductivity_dTPellets AuxKernel source file
	*
	*A method of calculating thermal conductivity
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "ThermalConductivity_dTPelletsAux.h"


template<>
InputParameters validParams<ThermalConductivity_dTPelletsAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_k_dT", "Set true to calculate Thermal conductivity of the pellets");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	params.addRequiredCoupledVar("burnup", "Coupled burnup");
	return params;
}

ThermalConductivity_dTPelletsAux::ThermalConductivity_dTPelletsAux(const std::string & name, InputParameters parameters)
  :AuxKernel(name, parameters),

	_model_k_dT(getParam<bool>("model_k_dT")),
  _temp(coupledValue("temp")),
	_burnup(coupledValue("burnup")),
	
	_porosity(getMaterialProperty<Real>("porosity"))
{
}

Real
ThermalConductivity_dTPelletsAux::computeValue()
{
	if (_model_k_dT)
	{
		const Real temp_ = _temp[_qp];
		const Real burnup_ = _burnup[_qp];
		const Real porosity_ = _porosity[_qp];

		//Thermal conductivity of UO2
		Real cond0, cond0_dT;
		Real kappa1d, kappa1d_dT;
		Real kappa1p, kappa1p_dT;
		Real kappa2p, kappa2p_dT;
		Real kappa4r, kappa4r_dT;

		computeChassieUnirradThCond(temp_, cond0, cond0_dT);
		computeDissFissionProdThCond(temp_, burnup_, kappa1d, kappa1d_dT);
		computePrecFissionProdThCond(temp_, burnup_, kappa1p, kappa1p_dT);
		computePoresThCond(temp_, porosity_, kappa2p, kappa2p_dT);
		computeRadDamageThCond(temp_, kappa4r, kappa4r_dT);

		const Real k_dT = cond0_dT * kappa1d * kappa1p * kappa2p * kappa4r + kappa1d_dT * cond0 * kappa1p * kappa2p * kappa4r + kappa1p_dT * cond0 * kappa1d * kappa2p * kappa4r + kappa2p_dT * cond0 * kappa1d * kappa1p * kappa4r +  kappa4r_dT * cond0 * kappa1d * kappa1p * kappa2p;

		return k_dT;
	}
	else
	{
		return 1; //default value
	}

}

//Compute Unirradiated Thermal Conductivity
void
ThermalConductivity_dTPelletsAux::computeChassieUnirradThCond(const Real temp, Real & cond0, Real & cond0_dT)
{
	const Real A = 0.030771;			//(m K)/W
	const Real B = 2.25e-4;				//m/W
	const Real C = 9.28e9;				//(W K)/m
	const Real D = 18295.09;			//K

	const Real temp2 = temp * temp;
	const Real temp4 = temp2 * temp2;
	
	cond0 = (1 /(A + B * temp)) + (C / temp2) * std::exp(-D / temp);
	cond0_dT = (C * std::exp(-D / temp) * (D - 2 * temp)) / (temp4) - B / (std::pow(A + B * temp, 2));
}

//Compute contribution to thermal conductivity due to dissolved fission products
void
ThermalConductivity_dTPelletsAux::computeDissFissionProdThCond(const Real temp, const Real burnup, Real & kappa1d, Real & kappa1d_dT)
{
	Real n = 0.0;
	Real m = 0.0;
	Real beta = burnup / 225;	//convert to atom percent

	if (beta <=0.1)
	{	
		kappa1d = 1.0;
		kappa1d_dT = 0.0;
	}
	else
	{
		beta = burnup / 225.;	//convert to atom percent
		m = 1.09 / std::pow(beta, 3.265);
		n = 0.0643 / std::sqrt(beta);
	
		kappa1d = (m + n * std::sqrt(temp)) * std::atan(1 / (m + n * std::sqrt(temp)));
		kappa1d_dT = (n / (2 * std::sqrt(temp))) * (std::atan(1 / (m + n * std::sqrt(temp))) - 1 / ((1 / std::pow(m + n * std::sqrt(temp), 2) + 1) * (m + n * std::sqrt(temp))));	
	}
}

//Compute contribution to thermal conductivity due to solid fission products
void
ThermalConductivity_dTPelletsAux::computePrecFissionProdThCond(const Real temp, const Real burnup, Real & kappa1p, Real & kappa1p_dT)
{
	Real beta = burnup / 225;  //convert to atom percent
	
	kappa1p = 1 + ((0.019 * beta) / (3 - 0.019 * beta)) * (1 / (1 + std::exp(-(temp - 1200.) / 100.)));
	kappa1p_dT = (-0.19 * beta*std::exp((temp + 1200.) / 100.)) / (19. * beta - 3000. * std::pow((std::exp(temp / 100.) + std::exp(12.)), 2));
}

//Computer contribution to thermal conductivity due to pores and fission gas bubbles
void
ThermalConductivity_dTPelletsAux::computePoresThCond(const Real temp, const Real porosity, Real & kappa2p, Real & kappa2p_dT)
{
	kappa2p = 1.0 - (2.05 - 5.0e-4 * temp) * porosity;
	kappa2p_dT = 0.0005 * porosity;
}

//Computer contribution to thermal conductivity due to radiation damage
void
ThermalConductivity_dTPelletsAux::computeRadDamageThCond(const Real temp, Real & kappa4r, Real & kappa4r_dT)
{
	kappa4r = 1 - 0.2 / (1 + std::exp((temp - 900.) / 80.));
	kappa4r_dT = (192.2 * std::exp(temp / 80.)) / std::pow((std::exp(temp / 80.) + std::exp(45. / 4.)), 2);
}


