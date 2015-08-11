/*VStrain AuxKernel source file
	*
	*Calculates the volumetric strain
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "VStrainAux.h"

template<>
InputParameters validParams<VStrainAux>()
{
  InputParameters params = validParams<AuxKernel>();
	params.addRequiredCoupledVar("densification_fraction", "Coupled densification factor");
	params.addRequiredParam<Real>("initial_porosity", "initial porosity");
  params.addRequiredParam<bool>("model_vstrain", "Set true to calculate volumetric strain");
	return params;
}

VStrainAux::VStrainAux(const InputParameters & parameters)
  :AuxKernel(parameters),
	_densificationF(coupledValue("densification_fraction")),
	_initial_porosity(getParam<Real>("initial_porosity")),
	_model_vstrain(getParam<bool>("model_vstrain")),

	_SFP(getMaterialProperty<Real>("solid_products")), //SFP swelling strain
	_GFP(getMaterialProperty<Real>("gaseous_products")) //GFP swelling strain
{
}

Real
VStrainAux::computeValue()
{
	if (_model_vstrain)
	{
		const Real densification_strain = (1 - _initial_porosity) / (1 - _initial_porosity * (1 - _densificationF[_qp])) - 1; //eq: 5.50 from prudil

		const Real vstrain = _SFP[_qp] + _GFP[_qp] + densification_strain;

		return vstrain;
	}
	else
	{
		return 1;
	}
}
