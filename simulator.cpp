#include <sbml/SBMLTypes.h>
#include <iostream>
#include <string>
#include <map>
#include <cmath>

LIBSBML_CPP_NAMESPACE_USE

using std::cout;
using std::endl;
using std::string;
using std::map;
using std::stod;
using std::runtime_error;

class System {
 private:
  const Model* model;
  map<string, double> global_params;

  double evaluateAst(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstName(
      const string& name,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstPlus(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstMinus(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstTimes(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstDivide(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;
  double evaluateAstPow(
      const ASTNode* node,
      const map<string, double>& species,
      const map<string, double>& local_params) const;

 public:
  explicit System(const Model* model);
  const Model* getModel() const;
  void evaluate(
      const map<string, double>& species,
      map<string, double>& species_dxdt) const;
};

System::System(const Model* model)
    : model {model} {
  auto num_global_params = model->getNumParameters();
  for (auto i = 0; i < num_global_params; i++) {
    auto p = this->model->getParameter(i);
    this->global_params[p->getId()] = p->getValue();
  }
}

const Model* System::getModel() const {
  return this->model;
}

void System::evaluate(
    const map<string, double>& species,
    map<string, double>& species_dxdt) const {
  // initialize dxdt
  for (auto it = species_dxdt.begin(); it != species_dxdt.end(); it++) {
    species_dxdt[it->first] = 0;
  }

  // reactions
  auto num_reactions = this->model->getNumReactions();
  for (auto i = 0; i < num_reactions; i++) {
    auto r = this->model->getReaction(i);
    auto kl = r->getKineticLaw();

    // local parameters
    map<string, double> local_params;
    auto num_local_params = kl->getNumParameters();
    for (auto j = 0; j < num_local_params; j++) {
      auto p = kl->getParameter(j);
      local_params[p->getId()] = p->getValue();
    }

    // evaluate AST node
    auto math = kl->getMath();
    auto result = this->evaluateAst(math, species, local_params);

    // update dxdt (reactants)
    auto num_reactants = r->getNumReactants();
    for (auto j = 0; j < num_reactants; j++) {
      auto sr = r->getReactant(j);
      auto sid = sr->getSpecies();
      species_dxdt[sid] -= result;
    }

    // update dxdt (products)
    auto num_products = r->getNumProducts();
    for (auto j = 0; j < num_products; j++) {
      auto sr = r->getProduct(j);
      auto sid = sr->getSpecies();
      species_dxdt[sid] += result;
    }
  }
}

double System::evaluateAst(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  switch (node->getType()) {
    case AST_NAME:
      return this->evaluateAstName(node->getName(), species, local_params);
    case AST_INTEGER:
      return static_cast<double>(node->getInteger());
    case AST_REAL:
      return node->getReal();
    case AST_PLUS:
      return this->evaluateAstPlus(node, species, local_params);
    case AST_MINUS:
      return this->evaluateAstMinus(node, species, local_params);
    case AST_TIMES:
      return this->evaluateAstTimes(node, species, local_params);
    case AST_DIVIDE:
      return this->evaluateAstDivide(node, species, local_params);
    case AST_POWER:
    case AST_FUNCTION_POWER:
      return this->evaluateAstPow(node, species, local_params);
    default:
      throw runtime_error("unknown node type");
  }
}

double System::evaluateAstName(
    const string& name,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  if (species.find(name) != species.end()) {
    return species.at(name);
  }

  if (local_params.find(name) != local_params.end()) {
    return local_params.at(name);
  }

  if (this->global_params.find(name) != this->global_params.end()) {
    return this->global_params.at(name);
  }

  throw runtime_error("unknown name: " + name);
}

double System::evaluateAstPlus(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  auto ret = 0.0;

  auto num_children = node->getNumChildren();
  for (auto i = 0; i < num_children; i++) {
    ret += this->evaluateAst(node->getChild(i), species, local_params);
  }

  return ret;
}

double System::evaluateAstMinus(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  auto ret = 0.0;

  auto num_children = node->getNumChildren();
  for (auto i = 0; i < num_children; i++) {
    ret -= this->evaluateAst(node->getChild(i), species, local_params);
  }

  return ret;
}

double System::evaluateAstTimes(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  auto ret = 1.0;

  auto num_children = node->getNumChildren();
  for (auto i = 0; i < num_children; i++) {
    ret *= this->evaluateAst(node->getChild(i), species, local_params);
  }

  return ret;
}

double System::evaluateAstDivide(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  auto ret = this->evaluateAst(node->getChild(0), species, local_params);

  auto num_children = node->getNumChildren();
  for (auto i = 1; i < num_children; i++) {
    auto denom = this->evaluateAst(node->getChild(i), species, local_params);
    if (denom == 0.0) {
      throw runtime_error("divide by zero");
    }
    ret /= denom;
  }

  return ret;
}

double System::evaluateAstPow(
    const ASTNode* node,
    const map<string, double>& species,
    const map<string, double>& local_params) const {
  auto num_children = node->getNumChildren();
  if (num_children != 2) {
    throw runtime_error("pow not binary");
  }

  auto left = this->evaluateAst(node->getLeftChild(), species, local_params);
  auto right = this->evaluateAst(node->getRightChild(), species, local_params);

  return pow(left, right);
}

class Stepper {
 public:
  void step(
      const System& system,
      const map<string, double>& species,
      map<string, double>& species_dx,
      double dt) const;
};

void Stepper::step(
    const System& system,
    const map<string, double>& species,
    map<string, double>& species_dx,
    double dt) const {
  map<string, double> species_dxdt;
  system.evaluate(species, species_dxdt);

  for (auto it = species_dxdt.begin(); it != species_dxdt.end(); it++) {
    species_dx[it->first] = it->second * dt;
  }
}

class Observer {
 public:
  void observe(double t, const map<string, double>& species) const;
};

void Observer::observe(double t, const map<string, double>& species) const {
  cout << t;
  for (auto it = species.begin(); it != species.end(); it++) {
    cout << " " << it->second;
  }
  cout << endl;
}

class Integrator {
 public:
  void integrate(
      const System& system,
      const Stepper& stepper,
      const Observer& observer,
      double duration,
      double dt) const;
};

void Integrator::integrate(
    const System& system,
    const Stepper& stepper,
    const Observer& observer,
    double duration,
    double dt) const {
  map<string, double> species;
  map<string, double> species_dx;

  // initialize species
  auto model = system.getModel();
  auto num_species = model->getNumSpecies();
  for (auto i = 0; i < num_species; i++) {
    auto s = model->getSpecies(i);
    species[s->getId()] = s->getInitialAmount();
  }

  for (auto t = 0.0; t <= duration; t += dt) {
    observer.observe(t, species);

    stepper.step(system, species, species_dx, dt);

    for (auto it = species_dx.begin(); it != species_dx.end(); it++) {
      species[it->first] += it->second;
    }
  }
}

void usage(char* cmd) {
  cout << "Usage: " << cmd << " <file> <duration> <dt>" << endl;
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
  if (argc < 4) {
    usage(argv[0]);
  }

  string filepath {argv[1]};
  auto duration = stod(argv[2]);
  auto dt = stod(argv[3]);

  SBMLReader reader;
  auto doc = reader.readSBMLFromFile(filepath);
  auto model = doc->getModel();

  System system(model);
  Stepper stepper;
  Observer observer;
  Integrator integrator;

  integrator.integrate(system, stepper, observer, duration, dt);

  return 0;
}
