#include <config.h>
#include "DWiener.h"
#include "DWieners.h"

#include <distribution/ScalarDist.h>

#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>

using std::vector;
using std::log;
using std::min;
using std::max;

#define BOUND(par) (*par[0])
#define TER(par) (*par[1])
#define BIAS(par) (*par[2])
#define DRIFT(par) (*par[3])
#define SD(par) (*par[3])

namespace wiener {

DWieners::DWieners(DWiener const *dist)
  : DWiener("dwieners", 5), _dist(dist)
{}

vector<double const *> DWieners::get_params(vector<double const *> const &par) const
{
  vector<double *> params(4);
  *params[0] = BOUND(par)*SD(par);
  *params[1] = TER(par);
  *params[2] = BIAS(par)*SD(par);
  *params[3] = DRIFT(par)*SD(par);

  vector<double const *> paramsconst(4);
  paramsconst[0] = params[0];
  paramsconst[1] = params[1];
  paramsconst[2] = params[2];
  paramsconst[3] = params[3];

  return paramsconst;
}

double DWieners::d(double x, PDFType type, 
         vector<double const *> const &par, bool give_log) const
{
  return _dist->d(x, type, get_params(par), give_log);
}

double DWieners::p_full(double q, std::vector<double const *> const &par, bool lower, 
  bool give_log) const
{
  return _dist->p_full(q, get_params(par), lower, give_log);
}

double DWieners::p(double q, std::vector<double const *> const &parameters, bool lower,
  bool give_log) const
{
  return _dist->p(q, get_params(parameters), lower, give_log);
}

double DWieners::q_full(double p, std::vector<double const *> const &parameters, bool lower,
  bool log_p) const
{
  return _dist->q_full(p, get_params(parameters), lower, log_p);
}

double DWieners::q(double p, std::vector<double const *> const &parameters, bool lower,
  bool log_p) const
{
  return _dist->q(p, get_params(parameters), lower, log_p);
}

double DWieners::r(std::vector<double const *> const &parameters, RNG *rng) const
{
  return _dist->r(get_params(parameters), rng);
}

}
