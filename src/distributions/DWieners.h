#ifndef DWIENERS_H_
#define DWIENERS_H_

#include "DWiener.h"

namespace jags {
namespace wiener {

/*
 * This class redefines the boundary separation parameter (alpha) and the drift
 * rate parameter (delta) and calls the apropriate d,p,q,r functions of the wiener
 * distribution functions with s=1.
 * Therefore only the d,p,q,r functions need to be overloaded, all other
 * functions can just be inherited.
 */
class DWieners : public DWiener
{
  DWiener const *_dist;

  public:
    DWieners(DWiener const *dist);

    bool checkParameterValue (std::vector<double const *> const &par) const;

    void get_params(std::vector<double const *> const &par, std::vector<double const *> &params) const;

    double d(double x, PDFType type,
      std::vector<double const *> const &parameters,
      bool give_log) const;

    double p_full(double q, std::vector<double const *> const &par, bool lower, 
      bool give_log) const;
    double p(double q, std::vector<double const *> const &parameters, bool lower,
      bool give_log) const;

    double q_full(double p, std::vector<double const *> const &parameters, bool lower,
      bool log_p) const;
    double q(double p, std::vector<double const *> const &parameters, bool lower,
      bool log_p) const;

    double r(std::vector<double const *> const &parameters, RNG *rng) const;
};

} //namespace wiener
} //namespace jags

#endif /* DWIENERS_H_ */
