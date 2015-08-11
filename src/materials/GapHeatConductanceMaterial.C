#include "GapHeatConductanceMaterial.h"

// Moose Includes
#include "PenetrationLocator.h"

// libMesh Includes
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<GapHeatConductanceMaterial>()
{
  MooseEnum orders("FIRST SECOND THIRD FOURTH", "FIRST");

  InputParameters params = validParams<Material>();
  params.addParam<std::string>("appended_property_name", "", "Name appended to material properties to make them unique");

  params.addRequiredCoupledVar("variable", "Temperature variable");
  params.addRequiredCoupledVar("k_pellets", "thermal conductivity of the pellets");
  params.addRequiredCoupledVar("k_sheath", "thermal conductivity of the sheath");
  // Node based
  params.addCoupledVar("gap_distance", "Distance across the gap");
  params.addCoupledVar("gap_temp", "Temperature on the other side of the gap");
  params.addParam<Real>("gap_conductivity", 1.0, "The thermal conductivity of the gap material");
  params.addParam<FunctionName>("gap_conductivity_function", "Thermal conductivity of the gap material as a function.  Multiplied by gap_conductivity.");
  params.addCoupledVar("gap_conductivity_function_variable", "Variable to be used in the gap_conductivity_function in place of time");


  // Quadrature based
  params.addParam<bool>("quadrature", false, "Whether or not to do quadrature point based gap heat transfer.  If this is true then gap_distance and gap_temp should NOT be provided (and will be ignored); however, paired_boundary and variable are then required.");
  params.addParam<BoundaryName>("paired_boundary", "The boundary to be penetrated");
  params.addParam<MooseEnum>("order", orders, "The finite element order");
  params.addParam<bool>("warnings", false, "Whether to output warning messages concerning nodes not being found");

  // Common
  params.addParam<Real>("min_gap", 1e-6, "A minimum gap size");
  params.addParam<Real>("max_gap", 1e6, "A maximum gap size");
  params.addParam<bool>("cylindrical_gap", false, "Use cylindrical wall heat flux");
  params.addParam<Real>("min_denom", 1e-6, "A minimum value for r*(r2/r1) term in cylindrical wall heat flux");
  params.addParam<Real>("max_denom", 1e6, "A maximum value for r*(r2/r1) term in cylindrical wall heat flux");

  params.addRequiredParam<Real>("emissivity_fuel", "The emissivity of the fuel surface");
  params.addRequiredParam<Real>("emissivity_sheath", "The emissivity of the cladding surface");

  params.addParam<bool>("use_displaced_mesh", true, "Whether or not this object should use the displaced mesh for computation.  Note that in the case this is true but no displacements are provided in the Mesh block the undisplaced mesh will still be used.");
	params.addRequiredParam<bool>("display_values", "display values");
  return params;
}

GapHeatConductanceMaterial::GapHeatConductanceMaterial(const InputParameters & parameters)
  :Material(parameters),
   _appended_property_name( getParam<std::string>("appended_property_name") ),
   _temp(coupledValue("variable")),
   _k_f(coupledValue("k_pellets")),
   _k_s(coupledValue("k_sheath")),
   _quadrature(getParam<bool>("quadrature")),
   _gap_temp(0),
   _gap_distance(88888),
   _radius(0),
   _r1(0),
   _r2(0),
   _has_info(false),
   _gap_distance_value(_quadrature ? _zero : coupledValue("gap_distance")),
   _gap_temp_value(_quadrature ? _zero : coupledValue("gap_temp")),
   _gap_conductance(declareProperty<Real>("gap_conductance"+_appended_property_name)),
   _gap_conductance_dT(declareProperty<Real>("gap_conductance"+_appended_property_name+"_dT")),
   _gap_conductivity(getParam<Real>("gap_conductivity")),
   _gap_conductivity_function(isParamValid("gap_conductivity_function") ? &getFunction("gap_conductivity_function") : NULL),
   _gap_conductivity_function_variable(isCoupled("gap_conductivity_function_variable") ? &coupledValue("gap_conductivity_function_variable") : NULL),
   _emissivity_f(getParam<Real>("emissivity_fuel")),
   _emissivity_s(getParam<Real>("emissivity_sheath")),
   _min_gap(getParam<Real>("min_gap")),
   _max_gap(getParam<Real>("max_gap")),
   _cylindrical_gap(getParam<bool>("cylindrical_gap")),
   _min_denom(getParam<Real>("min_denom")),
   _max_denom(getParam<Real>("max_denom")),
   _temp_var(_quadrature ? getVar("variable",0) : NULL),
   _penetration_locator(NULL),
   _serialized_solution(_quadrature ? &_temp_var->sys().currentSolution() : NULL),
   _dof_map(_quadrature ? &_temp_var->sys().dofMap() : NULL),
   _warnings(getParam<bool>("warnings")),
	 _display_values(getParam<bool>("display_values")),
	 _stefanboltzmann(5.6704e-8)
{
  if (_quadrature)
  {
    if (!parameters.isParamValid("paired_boundary"))
      mooseError(std::string("No 'paired_boundary' provided for ") + _name);
  }
  else
  {
    if (!isCoupled("gap_distance"))
      mooseError(std::string("No 'gap_distance' provided for ") + _name);

    if (!isCoupled("gap_temp"))
      mooseError(std::string("No 'gap_temp' provided for ") + _name);
  }


  if (_quadrature)
  {
    _penetration_locator = &_subproblem.geomSearchData().getQuadraturePenetrationLocator(parameters.get<BoundaryName>("paired_boundary"),
                                                                                         getParam<std::vector<BoundaryName> >("boundary")[0],
                                                                                         Utility::string_to_enum<Order>(parameters.get<MooseEnum>("order")));
  }
}


void
GapHeatConductanceMaterial::computeQpProperties()
{
  computeGapValues();
  computeQpConductance();
}


void
GapHeatConductanceMaterial::computeQpConductance()
{
  _gap_conductance[_qp] = h_gas_conduction() + h_solid_conduction() + h_radiation();
  _gap_conductance_dT[_qp] = dh_conduction() +dh_radiation();
}

//calculcats thermal conductivity of the pellet and sheath across the pellet - sheath gap
Real
GapHeatConductanceMaterial::h_solid_conduction()
{
	const Real a0 = 8.6e-6;
	const Real k_f(_k_f[_qp]);
	const Real k_s(_k_s[_qp]);
	const Real P_i(1e7); //need value
	const Real H = std::exp(26.034 - 0.026394 * _temp[_qp] + 4.3504e-5 * _temp[_qp] * _temp[_qp] * -2.5621e-8 * std::pow(_temp[_qp], 3));

	const Real h_solid_cond = ( (2 * k_f * k_s) / (k_f + k_s) ) * 1 / (a0 * H) * std::sqrt(P_i / _gap_distance); //eq: 5.5 from prudil
	return h_solid_cond;
}

//calculcats thermal conductivity of the gas in the sheath across the pellet - sheath gap
Real
GapHeatConductanceMaterial::h_gas_conduction()
{
	const Real k_g(2.5e-3); //thermal conductivity of the gas, needs to be updated after fission gas release is implemented
	const Real R_f(1e-6); //surface roughness of the fuel
	const Real R_s(5e-7); //surface roughness of the sheath
	const Real g = 5.20e-6; //jump distance of the sheath, = g_f + g_s, see p. 69

	const Real h_gas_cond = k_g / (1.5 * (R_f + R_s) + _gap_distance + g); //eq: 5.4 from prudil
	return h_gas_cond;
}

Real
GapHeatConductanceMaterial::dh_conduction()
{
  return 0;
}

//calculcats thermal conductivity of radiation across the pellet - sheath gap
Real
GapHeatConductanceMaterial::h_radiation()
{
  const Real numerator = (_temp[_qp]*_temp[_qp] + _gap_temp*_gap_temp) * (_temp[_qp] + _gap_temp);
	const Real denominator = 1 / _emissivity_f + 1 / _emissivity_s - 1;
  
return _stefanboltzmann * numerator / denominator; //eq: 5.18 from prudil
}

Real
GapHeatConductanceMaterial::dh_radiation()
{
  const Real numerator = 3 * _temp[_qp]*_temp[_qp] + _gap_temp * ( 2 * _temp[_qp] + _gap_temp );
	const Real denominator = 1 / _emissivity_f + 1 / _emissivity_s - 1;
  
	return _stefanboltzmann * numerator / denominator;
}

Real
GapHeatConductanceMaterial::gapLength(Real distance, Real min_gap, Real max_gap)
{
  Real gap_L = distance;

  if (gap_L > max_gap)
  {
    gap_L = max_gap;
  }

  gap_L = std::max(min_gap, gap_L);

  return gap_L;
}

Real
GapHeatConductanceMaterial::gapCyl(Real radius, Real r1, Real r2, Real min_denom, Real max_denom)
{
  Real denominator = radius*std::log(r2/r1);

  if (denominator > max_denom)
  {
    denominator = max_denom;
  }
  else if (denominator < min_denom)
  {
    denominator =  min_denom;
  }

  return denominator;
}


Real
GapHeatConductanceMaterial::gapK()
{
  Real gap_conductivity = _gap_conductivity;
  if (_gap_conductivity_function)
  {
    if (_gap_conductivity_function_variable)
    {
      gap_conductivity *= _gap_conductivity_function->value( (*_gap_conductivity_function_variable)[_qp], _q_point[_qp] );
    }
    else
    {
      gap_conductivity *= _gap_conductivity_function->value( _t, _q_point[_qp] );
    }
  }

  return gap_conductivity;
}

void
GapHeatConductanceMaterial::computeGapValues()
{
  if (!_quadrature)
  {
    _has_info = true;
    _gap_temp = _gap_temp_value[_qp];
    _gap_distance = _gap_distance_value[_qp];
  
		if (_display_values && _temp[_qp] > 1e3)
		{
		 	std::stringstream msg;
    	msg << "\tgap cond: " << _gap_conductance[_qp] << "\n"
					<< "\tgap temp: " << _gap_temp << "\n"
					<< "\tgap_distance: " << _gap_distance << "\n";
    	mooseWarning( msg.str() );
		}
}
  else
  {
    Node * qnode = _mesh.getQuadratureNode(_current_elem, _current_side, _qp);
    PenetrationInfo * pinfo = _penetration_locator->_penetration_info[qnode->id()];

    _gap_temp = 0.0;
    _gap_distance = 88888;
    _has_info = false;

    if (pinfo)
    {
      _gap_distance = pinfo->_distance;
      _has_info = true;

      Elem * slave_side = pinfo->_side;
      std::vector<std::vector<Real> > & slave_side_phi = pinfo->_side_phi;
      std::vector<unsigned int> slave_side_dof_indices;

      _dof_map->dof_indices(slave_side, slave_side_dof_indices, _temp_var->number());

      for (unsigned int i=0; i<slave_side_dof_indices.size(); ++i)
      {
        //The zero index is because we only have one point that the phis are evaluated at
        _gap_temp += slave_side_phi[i][0] * (*(*_serialized_solution))(slave_side_dof_indices[i]);
      }
    }
    else
    {
      if (_warnings)
        mooseWarning("No gap value information found for node " << qnode->id() << " on processor " << processor_id() << " at coordinate " << Point(*qnode));
    }
		if (_display_values && _temp[_qp] > 1e3)
		{
		 	std::stringstream msg;
    	msg << "\tgap cond: " << _gap_conductance[_qp] << "\n"
					<< "\tgap temp: " << _gap_temp << "\n"
					<< "\tgap_distance: " << _gap_distance << "\n";
    	mooseWarning( msg.str() );
		}
  }

  if (_cylindrical_gap)
  {
    if (_normals[_qp](0) > 0)
    {
      _r1 = _q_point[_qp](0);
      _r2 = _q_point[_qp](0) - _gap_distance; // note, _gap_distance is negative
      _radius = _r1;
    }
    else if (_normals[_qp](0) < 0)
    {
      _r1 = _q_point[_qp](0) + _gap_distance;
      _r2 = _q_point[_qp](0);
      _radius = _r2;
    }
    else
      mooseError( "Issue with cylindrical flux calc. normals. \n");
  }
}
