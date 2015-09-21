/*  This class contains functions for the wiener (diffusion) model */
#ifndef DWIENER_H_
#define DWIENER_H_

#include <distribution/ScalarDist.h>

#define WIENER_ERR 1e-10

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI  0.572364942924700087071713675677  /* log(sqrt(pi)) */
#endif

namespace jags {
namespace wiener {

/**
 * @short Diffusion Model distribution
 * <pre>
 * x ~ dwiener(v, a, w)
 * f(x | v, a, w) = ... if x slow
 *                = ... if x fast
 * </pre>
 * Parameters:
 * v - drift rate
 * a - boundary separation
 * w - relative starting point
 */
class DWiener : public ScalarDist 
{
  public:
    DWiener();
    DWiener(std::string const &name, unsigned int npar);

    /* 
     * logDensity, randomSample and typicalValue use the below defined
     * d,p,q,r functions
     */
    double logDensity(double x, PDFType type,
          std::vector<double const *> const &parameters,
          double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
       double const *lower, double const *upper,
       RNG *rng) const;
    double typicalValue(std::vector<double const *> const &parameters,
       double const *lower, double const *upper) const;
    /* 
     * Checks that:
     * a > 0
     * w in intervall [0,1] 
     * terr > 0
     */
    bool checkParameterValue(std::vector<double const *> const &parameters) const;

    // PDF function
    virtual double d(double x, PDFType type,
      std::vector<double const *> const &parameters,
      bool give_log) const;

    // CDF function
    virtual double p_full(double q, std::vector<double const *> const &par, bool lower, 
      bool give_log) const;
    virtual double p(double q, std::vector<double const *> const &parameters, bool lower,
      bool give_log) const;
    double F_lower(double q, double v, double a, double w) const;
    double Fl_lower(double q, double v, double a, double w, int K) const;
    double Fs_lower(double q, double v, double a, double w, int K) const;
    double Fs0_lower(double q, double a, double w, int K) const;
    double prob_upperbound(double v, double a, double w) const;
    double exp_pnorm(double a, double b) const;
    int K_small(double q, double v, double a, double w, double epsilon=WIENER_ERR) const;
    int K_large(double q, double v, double a, double w) const;

    // quantile function
    virtual double q_full(double p, std::vector<double const *> const &parameters, bool lower,
      bool log_p) const;
    virtual double q(double p, std::vector<double const *> const &parameters, bool lower,
      bool log_p) const;
    
    // random sampling function
    virtual double r(std::vector<double const *> const &parameters, RNG *rng) const;
    //double r_random_walk(std::vector<double const *> const &parameters, RNG *rng, double dt=0.0001) const;
    //double r_rejection_based(std::vector<double const *> const &parameters, RNG *rng) const;
    //double r_rejection_based_symmetric(std::vector<double> par, RNG *rng) const;
    double r_rejection_based_2(double a, double ter, double z, double v, RNG *rng) const;
};

} //namespace wiener
} //namespace jags

#endif /* DWIENER_H_ */
