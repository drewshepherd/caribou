#ifndef AVERAGEBURNUPAUX_H
#define AVERAGEBURNUPAUX_H

#include "AuxKernel.h"
#include "Function.h"
#include <math.h>

//Forward Declaration
class AverageBurnupAux;

template<>
InputParameters validParams<AuxKernel>();

class AverageBurnupAux : public AuxKernel
{
public:

  AverageBurnupAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

	Function * const _linear_power_time;
	Function * const _linear_power_burnup;

	const bool _model_wrt_time;
	const bool _model_wrt_burnup;

	const Real _density_initial;
private:
	const MaterialProperty<Real> & _ratio;
	const MaterialProperty<Real> & _pellet_radius;
	const MaterialProperty<Real> & _enrich;
 };

#endif //AVERAGEBURNUPAUX_H
