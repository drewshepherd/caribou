#ifndef CARIBOUAPP_H
#define CARIBOUAPP_H

#include "MooseApp.h"

class CaribouApp;

template<>
InputParameters validParams<CaribouApp>();

class CaribouApp : public MooseApp
{
public:
  CaribouApp(const std::string & name, InputParameters parameters);
  virtual ~CaribouApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* CARIBOUAPP_H */
