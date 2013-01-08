#ifndef DWIENERS_H_
#define DWIENERS_H_

#include "DWiener.h"

namespace wiener {

class DWieners : public DWiener
{
  DWiener const *_dist;

  public:
    DWieners(DWiener const *dist);

    std::vector<double const *> get_params(std::vector<double const *> const &par) const;

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

}

#endif /* DWIENERS_H_ */
