#ifndef DENSIFICATIONFRACTIONKERNEL_H
#define	DENSIFICATIONFRACTIONKERNEL_H

#include "Kernel.h"
#include "Function.h"
#include <math.h>

//Forward Declerations
class DensificationFractionKernel;

template<>
InputParameters validParams<DensificationFractionKernel>();

/**Temperature Dependent Material Properties of UO2**/
class DensificationFractionKernel : public Kernel
{
public:
  DensificationFractionKernel(const std::string & name,
      InputParameters parameters);

protected:
  virtual Real computeQpResidual();

	const bool _model_densification_fraction;
 
  const VariableValue  & _temp;
	const VariableValue  & _burnup_dt;
};
#endif //DENSIFICATIONFRACTIONKERNEL_H
