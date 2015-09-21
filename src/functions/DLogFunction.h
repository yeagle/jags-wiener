#ifndef DWIENER_DLOGFUNC_H_
#define DWIENER_DLOGFUNC_H_

#include "WFunction.h"

namespace jags {
namespace wiener {

class DLogFunction : public WFunction 
{
  public:
    DLogFunction(DWiener const *dist);

    double evaluate(std::vector<double const *> const &args) const;
};

} //namespace wiener
} //namespace jags

#endif /* DWIENER_DLOGFUNC_H_ */
