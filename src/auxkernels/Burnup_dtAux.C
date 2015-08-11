/*Burnup_dt AuxKernel source file
	*
	*the time derivative of BurnupAux
	*
	*written by Kyle Gamble and Drew Shepherd
*/

#include "Burnup_dtAux.h"
#include "Function.h"
#include <math.h>

template<>
InputParameters validParams<Burnup_dtAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_burnup_dt", "Set true to calculate timed derivative of burnup");
	return params;
}

Burnup_dtAux::Burnup_dtAux(const InputParameters & parameters)
  :AuxKernel(parameters),

	_model_burnup_dt(getParam<bool>("model_burnup_dt")),

  _density(getMaterialProperty<Real>("density")),
  _enrichment(getMaterialProperty<Real>("enrich")),
	_q_fission(getMaterialProperty<Real>("q_fission")),
	_ratio(getMaterialProperty<Real>("_ratio"))
{
}

Real
Burnup_dtAux::computeValue()
{
	if (_model_burnup_dt)
	{
		const Real MU235 = 235.0439; //mass of U-235
		const Real MU238 = 238.0508; //mass of U-238
		const Real MUO2 = 270.03; //mass of UO2
		const Real density_ = _density[_qp];
		const Real enrichment_ = _enrichment[_qp];
		const Real ratio_ = _ratio[_qp];
		const Real q_fission_ = _q_fission[_qp];

		//converts from density of uranium oxide to density of uranium
		const Real density_U = density_ * (MU235 * enrichment_ + MU238 * (1 - enrichment_)) / MUO2;

		//integration of q_fission over time (value of burnup during each time interval), [MWh/(kgU*s)]
		const Real value = (q_fission_ / (3.6e9 * ratio_ * density_U)); //eq: 5.26 from Prudil

		return value; 
	}
	else
	{
		return 1; //default value
	}
}


