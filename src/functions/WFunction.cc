#include <config.h>
#include "distributions/DWiener.h"
#include "WFunction.h"

using std::vector;
using std::string;

namespace jags {
namespace wiener {

WFunction::WFunction(string const &name, DWiener const *dist)
  :ScalarFunction(name, 5), _dist(dist)
{}

DWiener const *WFunction::dist() const
{
  return _dist;
}

bool WFunction::checkParameterValue(vector<double const *> const &args) const
{
  return checkArgs(args);
}

bool WFunction::checkArgs(vector<double const *> const &args) const
{
  vector<double const *> param(4);
  for (unsigned int i = 0; i < param.size(); ++i) {
    param[i] = args[i+1];
  }
  return _dist->checkParameterValue(param);
}

} //namespace wiener
} //namespace jags
