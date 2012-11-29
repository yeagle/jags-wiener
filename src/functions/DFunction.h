#ifndef DWIENER_DFUNC_H_
#define DWIENER_DFUNC_H_

#include <function/ScalarFunction.h>

namespace wiener {

class DFunction : public ScalarFunction 
{
  ScalarDist const *_dist;

  public:
    DFunction(std::string const &name, ScalarDist const *dist);

    ScalarDist const *dist() const;
    double evaluate(vector<double const *> const &args) const;
    bool checkParameterValue(vector<double const *> const &args) const;
    bool checkArgs(vector<double const *> const &args) const;
};

}

#endif /* DWIENER_DFUNC_H_ */
