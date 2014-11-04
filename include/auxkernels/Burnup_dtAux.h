#ifndef BURNUP_DTAUX_H
#define BURNUP_DTAUX_H

#include "AuxKernel.h"
#include "Function.h"
#include <math.h>

//Forward Declarations
class Burnup_dtAux;

template<>
InputParameters validParams<AuxKernel>();

class Burnup_dtAux : public AuxKernel
{
public:

  Burnup_dtAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

	const bool _model_burnup_dt;

private:
	const MaterialProperty<Real> & _density;
	const MaterialProperty<Real> & _enrichment;
	const MaterialProperty<Real> & _q_fission;
	const MaterialProperty<Real> & _ratio;
 };

#endif //BURNUP_DTAUX_H
