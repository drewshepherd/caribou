#ifndef THERMALCONDUCTIVITYSHEATHAUX_H
#define THERMALCONDUCTIVITYSHEATHAUX_H

#include "AuxKernel.h"


//Forward Declarations
class ThermalConductivitySheathAux;

template<>
InputParameters validParams<AuxKernel>();

class ThermalConductivitySheathAux : public AuxKernel
{
public:

  ThermalConductivitySheathAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

	const bool _model_k;
	const VariableValue & _temp;
private:

 };

#endif //THERMALCONDUCTIVITYSHEATHAUX_H
