#include <config.h>
#include "DFunction.h"

using std::vector;
using std::string;

namespace wiener {

DFunction::DFunction(ScalarDist const *dist)
  :ScalarFunction(dist->name, dist->npar() + 1), _dist(dist)
{}

ScalarDist const *DFunction::dist() const
{
  return _dist;
}

double DFunction::evaluate(vector<double const *> const &args) const
{
  double x = *args[0];
  vector<double const *> param(args.size() - 1);
  for (unsigned int i = 1; i < args.size(); ++i) {
    param[i-1] = args[i];
  }
  return dist()->d(x, PDF_FULL, param, false);
}

bool DFunction::checkParameterValue(vector<double const *> const &args) const
{
  if (dist()->discrete()) {
    double x = *args[0];
    if (x != static_cast<int>(x))
  return false;
  }
  return checkArgs(args);
}

bool DFunction::checkArgs(vector<double const *> const &args) const
{
  vector<double const *> param(_dist->npar());
  for (unsigned int i = 0; i < param.size(); ++i) {
    param[i] = args[i+1];
  }
  return _dist->checkParameterValue(param);
}

}
