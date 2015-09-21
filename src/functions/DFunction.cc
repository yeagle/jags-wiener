#include <config.h>
#include "distributions/DWiener.h"
#include "DFunction.h"

using std::vector;
using std::string;

namespace jags {
namespace wiener {

DFunction::DFunction(DWiener const *dist)
  :WFunction("dwiener", dist)
{}

double DFunction::evaluate(vector<double const *> const &args) const
{
  double x = *args[0];
  vector<double const *> param(args.size() - 1);
  for (unsigned int i = 1; i < args.size(); ++i) {
    param[i-1] = args[i];
  }
  return dist()->d(x, PDF_FULL, param, false);
}

} //namespace wiener
} //namespace jags    
