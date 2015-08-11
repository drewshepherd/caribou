/*ThermalConductivity_dTSheath AuxKernel source file
	*
	*A method of calculating thermal conductivity
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "ThermalConductivity_dTSheathAux.h"


template<>
InputParameters validParams<ThermalConductivity_dTSheathAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_k_dT", "Set true to calculate Thermal conductivity of the sheath");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	return params;
}

ThermalConductivity_dTSheathAux::ThermalConductivity_dTSheathAux(const InputParameters & parameters)
  :AuxKernel(parameters),

	_model_k_dT(getParam<bool>("model_k_dT")),
  _temp(coupledValue("temp"))
{
}

Real
ThermalConductivity_dTSheathAux::computeValue()
{
	if (_model_k_dT)
	{
		const Real temp2 = _temp[_qp] * _temp[_qp];

		const Real k_dT = 2.09e-2 - 2.9e-5 * _temp[_qp] + 23.07e-9 * temp2; //eq: 6.68 from Prudil

		return k_dT;
	}
	else
	{
		return 1; //default value
	}

}


