#ifndef THERMALCONDUCTIVITY_dTSHEATHAUX_H
#define THERMALCONDUCTIVITY_dTSHEATHAUX_H

#include "AuxKernel.h"


//Forward Declarations
class ThermalConductivity_dTSheathAux;

template<>
InputParameters validParams<AuxKernel>();

class ThermalConductivity_dTSheathAux : public AuxKernel
{
public:

  ThermalConductivity_dTSheathAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

	const bool _model_k_dT;
	const VariableValue & _temp;
private:

 };

#endif //THERMALCONDUCTIVITY_dTSHEATHAUX_H
