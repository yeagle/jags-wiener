#include <config.h>
#include "DLogFunction.h"

using std::vector;
using std::string;

namespace wiener {

DLogFunction::DLogFunction(ScalarDist const *dist)
  :ScalarFunction("dlogwiener", dist->npar() + 1), _dist(dist)
{}

ScalarDist const *DFunction::dist() const
{
  return _dist;
}

double DLogFunction::evaluate(vector<double const *> const &args) const
{
  double x = *args[0];
  vector<double const *> param(args.size() - 1);
  for (unsigned int i = 1; i < args.size(); ++i) {
    param[i-1] = args[i];
  }
  return dist()->d(x, PDF_FULL, param, true);
}

bool DLogFunction::checkParameterValue(vector<double const *> const &args) const
{
  if (dist()->discrete()) {
    double x = *args[0];
    if (x != static_cast<int>(x))
  return false;
  }
  return checkArgs(args);
}

bool DLogFunction::checkArgs(vector<double const *> const &args) const
{
  vector<double const *> param(_dist->npar());
  for (unsigned int i = 0; i < param.size(); ++i) {
    param[i] = args[i+1];
  }
  return _dist->checkParameterValue(param);
}

}
