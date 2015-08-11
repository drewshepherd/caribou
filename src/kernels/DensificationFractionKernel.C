/*DensificationFractionKernel source file
	*	
	*solves the densification fraction equation for the RHS
	*gives a percent from 0 - 60 as a DECIMAL VALUE (porosity is only 2 - 3 %, so the formula for the density change is (with percent converted to decimals is (1 - porosity) / ((1 - porosity) * (1 - densification_fraction))
	*
	*written by Drew Shepherd
*/

#include "DensificationFractionKernel.h"
#include "Function.h"
#include <math.h>

template<>
InputParameters validParams<DensificationFractionKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<bool>("model_densification_fraction", "Set to true to solve for densification fraction");
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
	params.addRequiredCoupledVar("burnup_dt", "Coupled burnup time derivative");
  return params;
}

DensificationFractionKernel::DensificationFractionKernel(const InputParameters & parameters) 
:Kernel(parameters),
	
	_model_densification_fraction(getParam<bool>("model_densification_fraction")),
  _temp(coupledValue("temp")),
	_burnup_dt(coupledValue("burnup_dt"))
{
}


//compute the densification factor that effects porosity, Fpower is the fraction of initial porosity which has been removed from the fuel
Real
DensificationFractionKernel::computeQpResidual()
{
	if (_model_densification_fraction)
	{
		const Real cd = 2.867e-2; // [kg/MWH]
		const Real bd = 8.67e-10; // [K^-3]
		const Real temp_ = _temp[_qp];
		const Real burnup_dt_ = _burnup_dt[_qp];
	
		Real temp_var = log(1 - _u[_qp] / 0.6) + bd * temp_ * temp_ * temp_;

			if (temp_var < 0)
			{
				temp_var = 0; //ensuring the log is not less than zero
			}

		const Real RHS = cd * (0.6 -_u[_qp]) * temp_var * burnup_dt_; //eg: 5.55 from Prudil, gives a unitless decimal
		return RHS * -_test[_i][_qp];
	}
	else
	{
		return 1; //default value
	}
}
