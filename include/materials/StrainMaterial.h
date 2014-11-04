#ifndef STRAINMATERIAL_H
#define STRAINMATERIAL_H

#include "StrainMaterial.h"
#include "SymmTensor.h"
#include "ColumnMajorMatrix.h"
#include "SolidMechanicsMaterial.h"
#include "SymmIsotropicElasticityTensor.h"
#include "VolumetricModel.h"
#include "SolidModel.h"

// Forward declarations
class StrainMaterial;

template<>
InputParameters validParams<StrainMaterial>();

class StrainMaterial : public SolidModel
{
public:
  StrainMaterial( const std::string & name, InputParameters parameters);
  virtual void modifyStrainIncrement();
	virtual Real computeGFP_dT(const Real density, const Real enrichment, const Real burnup, const Real burnup_dt, const Real temp);
	virtual void computeStress();
	virtual bool updateElasticityTensor(SymmElasticityTensor & tensor );
	virtual Real computeYoungsModulus(const Real temp);
private:
	const bool _model_thermal_expansion;
	const bool _model_youngs_modulus;
	const bool _display_values;

  const VariableValue  & _temp;
  const VariableValue  & _temp_old;
	const VariableValue  & _burnup;
	const VariableValue  & _burnup_old;
	const VariableValue  & _burnup_dt;
	const VariableValue  & _burnup_dt_old;
	const VariableValue  & _vstrain;
	const VariableValue  & _vstrain_old;

	Real _thermal_strain_increment;
	const MaterialProperty<Real> & _alpha1;
	const MaterialProperty<Real> & _alpha1_old;
	const MaterialProperty<Real> & _density;
	const MaterialProperty<Real> & _density_old;
	const MaterialProperty<Real> & _enrichment;
};

#endif //STRAINMATERIAL_H
