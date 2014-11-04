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

  ThermalConductivitySheathAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

	const bool _model_k;
	const VariableValue & _temp;
private:

 };

#endif //THERMALCONDUCTIVITYSHEATHAUX_H
