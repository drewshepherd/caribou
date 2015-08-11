/*Burnup AuxKernel source file
	*
	*A method of calculating burnup which differs from that of AverageBurnupAux
	*
	*written by Kyle Gamble and Drew Shepherd
*/
#include "BurnupAux.h"
#include "Function.h"
#include <math.h>

template<>
InputParameters validParams<BurnupAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_burnup", "Set true to calculate burnup");
	return params;
}

BurnupAux::BurnupAux(const InputParameters & parameters)
  :AuxKernel(parameters),

	_model_burnup(getParam<bool>("model_burnup")),

  _density(getMaterialProperty<Real>("density")),
  _enrichment(getMaterialProperty<Real>("enrich")),
	_q_fission(getMaterialProperty<Real>("q_fission")),
	_ratio(getMaterialProperty<Real>("_ratio"))
{
}

Real
BurnupAux::computeValue()
{
	if (_model_burnup)
	{
		const Real MU235 = 235.0439; //mass of U-235
		const Real MU238 = 238.0508; //mass of U-238
		const Real MUO2 = 270.03; //mass of UO2
		const Real density_ = _density[_qp];
		const Real enrichment_ = _enrichment[_qp] * 1e-2; //converts to a decimal
		const Real ratio_ = _ratio[_qp];
		const Real q_fission_ = _q_fission[_qp];

		//converts from density of uranium oxide to density of uranium
		const Real density_U = density_ * (MU235 * enrichment_ + MU238 * (1 - enrichment_)) / MUO2; 

		//integration of q_fission over time (value of burnup during each time interval), [MWh/(kgU)]
		const Real value = (q_fission_ / (3.6e3*1e6 * ratio_ * density_U)) * _dt; //eq: 5.26 from Prudil

		return _u_old[_qp] + value; //sum of burnup
	}
	else
	{
		return 1; //default value
	}

}


