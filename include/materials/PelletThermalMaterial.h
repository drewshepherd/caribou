#ifndef PELLETTHERMALMATERIAL_H
#define	PELLETTHERMALMATERIAL_H

#include "Material.h"

//Forward Declerations
class PelletThermalMaterial;

template<>
InputParameters validParams<PelletThermalMaterial>();

/**Temperature Dependent Material Properties of UO2**/
class PelletThermalMaterial : public Material
{
public:
  PelletThermalMaterial(const InputParameters & parameters);

  virtual void initStatefulProperties(unsigned n_points);

protected:
  virtual void computeProperties();
  virtual Real computeSpecificHeat(const Real temp);
  virtual Real computeAlpha(const Real temp);
  virtual Real computeSFP(const Real burnup);
  virtual Real computeGFP(const Real density, const Real enrichment, const Real burnup, const Real burnup_dt, const Real temp);

	Real computeTheoDensity(const Real alpha, const Real temp);
  Real computeSFPFactor(const Real burnup);

	const bool _model_thermal_conductivity;
	const bool _model_specific_heat;
	const bool _model_porosity;
	const bool _model_alpha;
	const bool _model_SFP;
	const bool _model_GFP;
	const bool _display_values;

  const VariableValue  & _temp;
	const VariableValue  & _burnup;
	const VariableValue  & _burnup_dt;
	const VariableValue  & _densificationF;
  const VariableValue    & _k;
  const VariableValue    & _k_dT;

	const Real _initial_porosity;

private:
	const MaterialProperty<Real> & _density;
	const MaterialProperty<Real> & _enrichment;
  MaterialProperty<Real> & _thermal_conductivity;
  MaterialProperty<Real> & _thermal_conductivity_dT;
  MaterialProperty<Real> & _specific_heat;
	MaterialProperty<Real> & _alpha1;
	MaterialProperty<Real> & _alpha1_old;
	MaterialProperty<Real> & _SFP;
	MaterialProperty<Real> & _GFP;
	MaterialProperty<Real> & _porosity;
};

#endif //PELLETTHERMALMATERIAL_H
