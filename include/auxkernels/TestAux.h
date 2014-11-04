#ifndef TESTAUX_H
#define TESTAUX_H

#include "AuxKernel.h"

//Forward Declarations
class TestAux;

template<>
InputParameters validParams<AuxKernel>();

class TestAux : public AuxKernel
{
public:

  TestAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
	const	bool _model_density;
	const bool _model_q_fission;
	const bool _model_q_fission_old;
	const bool _model_ratio;
	const bool _model_pellet_radius;
	const bool _model_porosity;
	const bool _model_enrichment;
	const bool _model_SFP;
	const bool _model_GFP;
private:
	const MaterialProperty<Real> & _density;
	const MaterialProperty<Real> & _q_fission;
  const MaterialProperty<Real> & _q_fission_old;
  const MaterialProperty<Real> & _ratio;
	const MaterialProperty<Real> & _pellet_radius;
	const MaterialProperty<Real> & _porosity;
  const MaterialProperty<Real> & _enrichment;
  const MaterialProperty<Real> & _SFP;
  const MaterialProperty<Real> & _GFP;
 };

#endif //TESTAUX_H
