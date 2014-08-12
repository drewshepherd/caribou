#include "CaribouApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<CaribouApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

CaribouApp::CaribouApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

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
}

void
CaribouApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
