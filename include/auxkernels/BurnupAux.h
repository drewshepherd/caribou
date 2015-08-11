#ifndef BURNUPAUX_H
#define BURNUPAUX_H

#include "AuxKernel.h"
#include "Function.h"
#include <math.h>

//Forward Declarations
class BurnupAux;

template<>
InputParameters validParams<AuxKernel>();

class BurnupAux : public AuxKernel
{
public:

  BurnupAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

	const bool _model_burnup;
private:
	const MaterialProperty<Real> & _density;
	const MaterialProperty<Real> & _enrichment;
	const MaterialProperty<Real> & _q_fission;
	const MaterialProperty<Real> & _ratio;
 };

#endif //BURNUPAUX_H
