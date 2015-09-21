#ifndef DWIENER_DFUNC_H_
#define DWIENER_DFUNC_H_

#include "WFunction.h"

namespace jags {
namespace wiener {

class DFunction : public WFunction 
{
  public:
    DFunction(DWiener const *dist);

    double evaluate(std::vector<double const *> const &args) const;
};

} //namespace wiener
} //namespace jags

#endif /* DWIENER_DFUNC_H_ */
