#ifndef SHEATHMECHANICALMATERIAL_H
#define SHEATHMECHANICALMATERIAL_H

#include "SolidModel.h"

// Forward declarations
class SheathMechanicalMaterial;

template<>
InputParameters validParams<SheathMechanicalMaterial>();

class SheathMechanicalMaterial : public SolidModel
{
public:
  SheathMechanicalMaterial(const InputParameters & parameters);

protected:
  const VariableValue  & _temp;

	const Real _relative_tolerance;
  const Real _absolute_tolerance;
  const unsigned int _max_its;
  const bool _output_iteration_info;

  const bool _model_thermal_expansion;
	const bool _model_youngs_modulus;
	const bool _model_diffusional_creep;
	const bool _display_values;

	MaterialProperty<SymmTensor> & _creep_strain;
  MaterialProperty<SymmTensor> & _creep_strain_old;
  MaterialProperty<Real> & _primary_creep_strain;
  MaterialProperty<Real> & _primary_creep_strain_old;
	
	virtual Real computeAxialThermEx(const Real temp);
	virtual Real computeRadialThermEx(const Real temp);
	virtual Real computeYoungsModulus(const Real temp);
  virtual void modifyStrainIncrement();
	virtual void computeStress();
	virtual bool updateElasticityTensor(SymmElasticityTensor & tensor );
	virtual double EDOTGB(double T, double sigma_a);
	virtual double d_EDOTGB_d_sigma(double T, double sigma_a);
};

#endif //SHEATHMECHANICALMATERIAL_H
