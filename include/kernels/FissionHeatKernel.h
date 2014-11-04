#ifndef FISSIONHEATKERNEL_H
#define FISSIONHEATKERNEL_H

#include "Kernel.h"
#include "Function.h"
#include <math.h>

//Forward Declarations
class FissionHeatKernel;

template<>
InputParameters validParams<FissionHeatKernel>();

class FissionHeatKernel : public Kernel
{
public:

  FissionHeatKernel(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();
	
private:
	const MaterialProperty<Real> & _q_fission;
 };

#endif //FISSIONHEATKERNEL_H
