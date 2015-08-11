#ifndef VSTRAINAUX_H
#define VSTRAINAUX_H

#include "AuxKernel.h"

//Forward Declaration
class VStrainAux;

template<>
InputParameters validParams<AuxKernel>();

class VStrainAux : public AuxKernel
{
public:

  VStrainAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

	const VariableValue  & _densificationF;

	const Real _initial_porosity;
	bool _model_vstrain;
private:
	const MaterialProperty<Real> & _SFP;
	const MaterialProperty<Real> & _GFP;
 };

#endif //VSTRAINAUX_H
