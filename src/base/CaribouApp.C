#include "CaribouApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

//AuxKernels
#include "AverageBurnupAux.h"
#include "BurnupAux.h"
#include "Burnup_dtAux.h"
#include "VStrainAux.h"
#include "ThermalConductivityPelletsAux.h"
#include "ThermalConductivity_dTPelletsAux.h"
#include "ThermalConductivitySheathAux.h"
#include "ThermalConductivity_dTSheathAux.h"
#include "TestAux.h"

//Kernels
#include "FissionHeatKernel.h"
#include "DensificationFractionKernel.h"

//Material
#include "FissionHeatMaterial.h"
#include "PelletThermalMaterial.h"
#include "SheathThermalMaterial.h"
#include "SheathMechanicalMaterial.h"
#include "StrainMaterial.h"
#include "GapHeatConductanceMaterial.h"

template<>
InputParameters validParams<CaribouApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

CaribouApp::CaribouApp(InputParameters parameters) :
    MooseApp(parameters)
{

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  CaribouApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  CaribouApp::associateSyntax(_syntax, _action_factory);
}

CaribouApp::~CaribouApp()
{
}

void
CaribouApp::registerApps()
{
  registerApp(CaribouApp);
}

void
CaribouApp::registerObjects(Factory & factory)
{
	//AuxKernels
	registerAux(AverageBurnupAux);
	registerAux(BurnupAux);
	registerAux(Burnup_dtAux);
	registerAux(VStrainAux);
	registerAux(ThermalConductivityPelletsAux);
	registerAux(ThermalConductivity_dTPelletsAux);
	registerAux(ThermalConductivitySheathAux);
	registerAux(ThermalConductivity_dTSheathAux);
	registerAux(TestAux);

	//Kernels
	registerKernel(FissionHeatKernel);
	registerKernel(DensificationFractionKernel);

	//Material
	registerMaterial(FissionHeatMaterial);
	registerMaterial(PelletThermalMaterial);
	registerMaterial(SheathThermalMaterial);
	registerMaterial(SheathMechanicalMaterial);
	registerMaterial(StrainMaterial);
	registerMaterial(GapHeatConductanceMaterial);
}

void
CaribouApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
