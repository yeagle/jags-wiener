#ifndef DWIENER_WFUNC_H_
#define DWIENER_WFUNC_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace wiener {

class WFunction : public ScalarFunction 
{
  DWiener const *_dist;

  public:
    WFunction(std::string const &name, DWiener const *dist);

    DWiener const *dist() const;
    virtual double evaluate(std::vector<double const *> const &args) const=0;
    bool checkParameterValue(std::vector<double const *> const &args) const;
    bool checkArgs(std::vector<double const *> const &args) const;
};

} //namespace wiener
} //namespace jags

#endif /* DWIENER_WFUNC_H_ */
