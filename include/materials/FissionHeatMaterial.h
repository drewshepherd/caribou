#ifndef FISSIONHEATMATERIAL_H
#define FISSIONHEATMATERIAL_H

#include "Material.h"
#include <math.h>
#include "Function.h"

//using namespace std;

//Forward declaration
class FissionHeatMaterial;

template <>
InputParameters validParams<FissionHeatMaterial>();

class FissionHeatMaterial : public Material
{
public:
	//Constructor - obtains the values of burnup, enrichment and pellet radius from the input file
	FissionHeatMaterial(const InputParameters & parameters);

protected:
	double interp2(const double _cons_array[][10]);
	virtual Real bessel0(double x);
	virtual Real bessel1(double x);
	Real Q_fission();
  virtual void computeQpProperties();
	virtual void initQpStatefulProperties();
	
	Function * const _linear_power;

	const bool _model_Qfission;

  const VariableValue & _burnup_avg;
	const Real _enrichment_property;
	const Real _pellet_radius_property;

	const Real _ratioProperty;
	const Real _initial_density;
	const Real _initial_qfission;
	const Real _area_fuel;

	const bool _is_3D;
	const bool _model_plate_fuel;
private:
	const MaterialProperty<Real> & _density;
	MaterialProperty<Real> & _q_fission;
	MaterialProperty<Real> & _q_fission_old;
  MaterialProperty<Real> & _ratio;
  MaterialProperty<Real> & _pellet_radius;
  MaterialProperty<Real> & _enrichment;
};
#endif //FISSIONHEATMATERIAL_H
