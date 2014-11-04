#ifndef GAPHEATCONDUCTANCEMATERIAL_H
#define GAPHEATCONDUCTANCEMATERIAL_H

#include "Material.h"

/**
 * Generic gap heat transfer model, with h_gap =  h_conduction + h_contact + h_radiation
 */
class GapHeatConductanceMaterial :
  public Material
{
public:

  GapHeatConductanceMaterial(const std::string & name, InputParameters parameters);

  virtual ~GapHeatConductanceMaterial(){}

  static Real gapLength(Real distance, Real min_gap, Real max_gap);

  static Real gapCyl( Real radius, Real r1, Real r2, Real min_denom, Real max_denom);

protected:

  virtual void computeQpProperties();

  /**
   * Override this to compute the conductance at _qp
   */
  virtual void computeQpConductance();

  virtual Real h_gas_conduction();
  virtual Real h_solid_conduction();
  virtual Real h_radiation();
  virtual Real dh_conduction();
  virtual Real dh_radiation();
  virtual Real gapK();

  virtual void computeGapValues();

  const std::string _appended_property_name;

  const VariableValue & _temp;
  const VariableValue & _k_f;
  const VariableValue & _k_s;

  bool _quadrature;

  Real _gap_temp;
  Real _gap_distance;
  Real _radius;
  Real _r1;
  Real _r2;

  bool _has_info;

  const VariableValue & _gap_distance_value;
  const VariableValue & _gap_temp_value;
  MaterialProperty<Real> & _gap_conductance;
  MaterialProperty<Real> & _gap_conductance_dT;

  const Real _gap_conductivity;
  Function * const _gap_conductivity_function;
  const VariableValue * _gap_conductivity_function_variable;

  const Real _emissivity_f;
  const Real _emissivity_s;

  Real _min_gap;
  Real _max_gap;
  const bool _cylindrical_gap;
  Real _min_denom;
  Real _max_denom;

  MooseVariable * _temp_var;
  PenetrationLocator * _penetration_locator;
  const NumericVector<Number> * * _serialized_solution;
  DofMap * _dof_map;
  const bool _warnings;
  const bool _display_values;
  const Real _stefanboltzmann;
};

template<>
InputParameters validParams<GapHeatConductanceMaterial>();

#endif //GAPHEATCONDUCTANCEMATERIAL_H
