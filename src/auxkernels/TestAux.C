/*TestAux AuxKernel source file
	*
	*Displays the values that cannot be seen directly in Paraview. Used for testing/troubleshooting purposes
	*
	*written by  by Drew Shepherd
*/

#include "TestAux.h"

template<>
InputParameters validParams<TestAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<bool>("model_density", "Set true to calculate density");
  params.addRequiredParam<bool>("model_q_fission", "Set true to calculate q_fission");
  params.addRequiredParam<bool>("model_q_fission_old", "Set true to calculate q_fission_old");
  params.addRequiredParam<bool>("model_ratio", "Set true to calculate ratio");
  params.addRequiredParam<bool>("model_pellet_radius", "Set true to calculate pellet_radius");
  params.addRequiredParam<bool>("model_porosity", "Set true to calculate porosity");
  params.addRequiredParam<bool>("model_enrichment", "Set true to calculate enrichment");
  params.addRequiredParam<bool>("model_SFP", "Set true to calculate solid fission product strain");
  params.addRequiredParam<bool>("model_GFP", "Set true to calculate gaseous fission product strain");

	return params;
}

TestAux::TestAux(const std::string & name, InputParameters parameters)
  :AuxKernel(name, parameters),

	_model_density(getParam<bool>("model_density")),
	_model_q_fission(getParam<bool>("model_q_fission")),
	_model_q_fission_old(getParam<bool>("model_q_fission_old")),
	_model_ratio(getParam<bool>("model_ratio")),
	_model_pellet_radius(getParam<bool>("model_pellet_radius")),
	_model_porosity(getParam<bool>("model_porosity")),
	_model_enrichment(getParam<bool>("model_enrichment")),
	_model_SFP(getParam<bool>("model_SFP")),
	_model_GFP(getParam<bool>("model_GFP")),

	_density(getMaterialProperty<Real>("density")),
	_q_fission(getMaterialProperty<Real>("q_fission")),
	_q_fission_old(getMaterialPropertyOld<Real>("q_fission")),
	_ratio(getMaterialProperty<Real>("_ratio")),
	_pellet_radius(getMaterialProperty<Real>("pellet_rad")),
	_porosity(getMaterialProperty<Real>("porosity")),
	_enrichment(getMaterialProperty<Real>("enrich")),
	_SFP(getMaterialProperty<Real>("solid_products")),
	_GFP(getMaterialProperty<Real>("gaseous_products"))
{}

Real
TestAux::computeValue()
{
	if (_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _density[_qp];
	}
	else if (!_model_density && _model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _q_fission[_qp];
	}
	else if (!_model_density && !_model_q_fission && _model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _q_fission_old[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && _model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _ratio[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && _model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _pellet_radius[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && _model_porosity && !_model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _porosity[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && _model_enrichment && !_model_SFP && !_model_GFP)
	{
		return _enrichment[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && _model_SFP && !_model_GFP)
	{
		return _SFP[_qp];
	}
	else if (!_model_density && !_model_q_fission && !_model_q_fission_old && !_model_ratio && !_model_pellet_radius && !_model_porosity && !_model_enrichment && !_model_SFP && _model_GFP)
	{
		return _GFP[_qp];
	}
	else
	{
		mooseError("There must be one true parameter and the other parameters all false for test_aux to run");
		return 1;
	}
}
