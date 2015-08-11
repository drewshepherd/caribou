/*FissionHeat Kernel source file
	*
	*This is one term of the heat equation. The other terms can be found in the input file.
	*multiplies q_fission by a test function, this allows us to use q_fission as a heat source term coupled to temp
	*
	*written by Kyle Gamble and by Drew Shepherd
*/

#include "FissionHeatKernel.h"

template<>
InputParameters validParams<FissionHeatKernel>()
{
  InputParameters params = validParams<Kernel>();

  return params;
}

FissionHeatKernel::FissionHeatKernel(const InputParameters & parameters)
  :Kernel(parameters),
	_q_fission(getMaterialProperty<Real>("q_fission"))
{
}

Real
FissionHeatKernel::computeQpResidual()
{
  return -_test[_i][_qp] * _q_fission[_qp]; //eq: 5.1 from Prudil
}
