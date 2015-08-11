/*ThermalConductivitySheath AuxKernel source file
	*
	*A method of calculating thermal conductivity
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "ThermalConductivitySheathAux.h"


template<>
InputParameters validParams<ThermalConductivitySheathAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_k", "Set true to calculate Thermal conductivity of the sheath");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	return params;
}

ThermalConductivitySheathAux::ThermalConductivitySheathAux(const InputParameters & parameters)
  :AuxKernel(parameters),

	_model_k(getParam<bool>("model_k")),
  _temp(coupledValue("temp"))
{
}

Real
ThermalConductivitySheathAux::computeValue()
{
	if (_model_k)
	{
		const Real temp2 = _temp[_qp] * _temp[_qp];
		const Real temp3 = temp2 * _temp[_qp];

		const Real k = 7.51 + 2.09e-2 * _temp[_qp] -1.45e-5 * temp2 + 7.67e-9 * temp3; //eq: 6.68 from Prudil
		
		return k;
	}
	else
	{
		return 1; //default value
	}

}


