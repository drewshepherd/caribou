#ifndef SHEATHTHERMALMATERIAL_H
#define SHEATHTHERMALMATERIAL_H

#include "Material.h"

//Forward Declarations
class SheathThermalMaterial;

template<>
InputParameters validParams<SheathThermalMaterial>();

/**
 * Temperature  dependent thermal properties of zirconium alloy
 */

class SheathThermalMaterial : public Material
{
public:
  SheathThermalMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();
	virtual Real computeSpecificHeat(const Real temp);

	const bool _model_thermal_conductivity;
	const bool _model_specific_heat;

  const VariableValue    & _temp;
  const VariableGradient & _grad_temp;
  const VariableValue    & _k;
  const VariableValue    & _k_dT;
private:

  MaterialProperty<Real>     & _thermal_conductivity;
  MaterialProperty<Real>     & _thermal_conductivity_dT;
  MaterialProperty<Real>     & _specific_heat;
};

#endif //SHEATHTHERMALMATERIAL_H
