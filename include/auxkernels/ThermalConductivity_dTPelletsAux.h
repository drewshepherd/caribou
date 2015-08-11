#ifndef THERMALCONDUCTIVITY_dTPELLETSAUX_H
#define THERMALCONDUCTIVITY_dTPELLETSAUX_H

#include "AuxKernel.h"


//Forward Declarations
class ThermalConductivity_dTPelletsAux;

template<>
InputParameters validParams<AuxKernel>();

class ThermalConductivity_dTPelletsAux : public AuxKernel
{
public:

  ThermalConductivity_dTPelletsAux(const InputParameters & parameters);
protected:
  virtual Real computeValue();
  virtual void computeChassieUnirradThCond(const Real temp, Real & cond0, Real & cond0_dT);
	virtual void computeDissFissionProdThCond(const Real temp, const Real burnup, Real & kappa1d, Real & kappa1d_dT);
	virtual void computePrecFissionProdThCond(const Real temp, const Real burnup, Real & kappa1p, Real & kappa1p_dT);
	virtual void computePoresThCond(const Real temp, const Real porosity, Real & kappa2p, Real & kappa2p_dT);
	virtual void computeRadDamageThCond(const Real temp, Real & kappa4r, Real & kappa4r_dT);

	const bool _model_k_dT;
	const VariableValue & _temp;
	const VariableValue & _burnup;
private:
	const MaterialProperty<Real> & _porosity;
 };

#endif //THERMALCONDUCTIVITY_dTPELLETSAUX_H
