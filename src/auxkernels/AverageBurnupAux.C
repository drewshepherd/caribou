/*Burnup AuxKernel source file
	*
	*gives the sum of burnup over all space, averaged over the volume
	*
	*written by Kyle Gamble and by Drew Shepherd
*/

#include "AverageBurnupAux.h"
#include "Function.h"

template<>
InputParameters validParams<AverageBurnupAux>()
{
  InputParameters params = validParams<AuxKernel>();

	params.addParam<FunctionName>("linear_power_time", "The linear element power function (W/m) as a function of time");
	params.addParam<FunctionName>("linear_power_burnup", "The linear element power function (W/m) as a function of burnup");
  params.addRequiredParam<bool>("model_wrt_time", "Set true to calculate average burnup as a function of time");
  params.addRequiredParam<bool>("model_wrt_burnup", "Set true to calculate average burnup as a function of burnup");
	params.addRequiredParam<Real>("initial_density", "The initial pellet density (kg/m^3)");
	return params;
}

AverageBurnupAux::AverageBurnupAux(const std::string & name, InputParameters parameters)
  :AuxKernel(name, parameters),

	_linear_power_time(&getFunction("linear_power_time")), //in units of [W/m]
	_linear_power_burnup(&getFunction("linear_power_burnup")), //in units of [W/m]
	_model_wrt_time(getParam<bool>("model_wrt_time")),
	_model_wrt_burnup(getParam<bool>("model_wrt_burnup")),
	_density_initial(getParam<Real>("initial_density")),

	_ratio(getMaterialProperty<Real>("_ratio")),
	_pellet_radius(getMaterialProperty<Real>("pellet_rad")),
	_enrich(getMaterialProperty<Real>("enrich"))
{
}

Real
AverageBurnupAux::computeValue()
{
	const double MUO2 = 270.03; //mass of urandium oxide
  const double pi = 3.141592654;
	const double MU235 = 235.0439; //mass of U-235
	const double MU238 = 238.0508; //mass of U-238
	const Real enrichment_ = _enrich[_qp] * 1e-2; //convert enrich to decimal
	const Real pellet_radius_ = _pellet_radius[_qp];
	const Real ratio_ = _ratio[_qp];

	const Real MUO2_L = _density_initial * pi * pellet_radius_ * pellet_radius_; //mass of urandium oxide per unit length
	const Real MU_L = (MU235 * enrichment_ + MU238 * (1 - enrichment_)) * MUO2_L / MUO2; //mass of uranium per unit length
	
	Real linear_power(0.);
	if (_model_wrt_time)
	{
		linear_power = _linear_power_time->value(_t, _q_point[_qp]);
	}
	else if (_model_wrt_burnup)
	{
		linear_power = _linear_power_burnup->value(_t, _q_point[_qp]);
	}
	else
	{
		mooseError("must specify either a time or burnup dependent linear power function");
		return 1; //default value
	}
	
	const Real value = linear_power/ (3.6e3 * 1e6 * MU_L * ratio_) * _dt; //eq: 5.27 from Prudil
	return _u_old[_qp] + value;
}
