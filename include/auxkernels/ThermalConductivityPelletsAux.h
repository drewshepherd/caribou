#ifndef THERMALCONDUCTIVITYPELLETSAUX_H
#define THERMALCONDUCTIVITYPELLETSAUX_H

#include "AuxKernel.h"


//Forward Declarations
class ThermalConductivityPelletsAux;

template<>
InputParameters validParams<AuxKernel>();

class ThermalConductivityPelletsAux : public AuxKernel
{
public:

  ThermalConductivityPelletsAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  virtual void computeChassieUnirradThCond(const Real temp, Real & cond0);
	virtual void computeDissFissionProdThCond(const Real temp, const Real burnup, Real & kappa1d);
	virtual void computePrecFissionProdThCond(const Real temp, const Real burnup, Real & kappa1p);
	virtual void computePoresThCond(const Real temp, const Real porosity, Real & kappa2p);
	virtual void computeRadDamageThCond(const Real temp, Real & kappa4r);

	const bool _model_k;
	const VariableValue & _temp;
	const VariableValue & _burnup;
private:
	const MaterialProperty<Real> & _porosity;
 };

#endif //THERMALCONDUCTIVITYPELLETSAUX_H
