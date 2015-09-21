#include <config.h>
#include "DWiener.h"
#include "DWieners.h"

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
#define SD(par) (*par[4])

/*
 * TODO : This file contains really bad code, that slows the wieners node extremely!
 * Contributions to improve this code are encouraged and very welcome.
 *
 * I just don't know how to alter the ' vector<double const *> const '
 * parameter vector in an efficient way...
 *
 * There's also a lot of redudnant code to get rid of...
 */

namespace jags {
namespace wiener {

DWieners::DWieners(DWiener const *dist)
  : DWiener("dwieners", 5), _dist(dist)
{}

bool DWieners::checkParameterValue (vector<double const *> const &par) const
{
    return (BOUND(par)*SD(par)>0 && 
            BIAS(par)<1 && BIAS(par)>0 &&
            TER(par) > 0 && SD(par) > 0);
}

void DWieners::get_params(vector<double const *> const &par, vector<double const *> &params) const
{
  // TODO : replace redundant code and use this function
}

double DWieners::d(double x, PDFType type, 
         vector<double const *> const &par, bool give_log) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;

  return _dist->d(x, type, params, give_log);
}

double DWieners::p_full(double q, vector<double const *> const &par, bool lower, 
  bool give_log) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;
  return _dist->p_full(q, params, lower, give_log);
}

double DWieners::p(double q, vector<double const *> const &par, bool lower,
  bool give_log) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;
  return _dist->p(q, params, lower, give_log);
}

double DWieners::q_full(double p, vector<double const *> const &par, bool lower,
  bool log_p) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;
  return _dist->q_full(p, params, lower, log_p);
}

double DWieners::q(double p, vector<double const *> const &par, bool lower,
  bool log_p) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;
  return _dist->q(p, params, lower, log_p);
}

double DWieners::r(vector<double const *> const &par, RNG *rng) const
{
  double bound,ter,bias,drift;
  bound = BOUND(par)/SD(par);
  ter = TER(par);
  bias = BIAS(par);
  drift = DRIFT(par)/SD(par);

  double *boundp = &bound; 
  double *terp = &ter; 
  double *biasp = &bias; 
  double *driftp = &drift; 

  vector<double const *> params(4);
  params[0] = boundp;
  params[1] = terp;
  params[2] = biasp;
  params[3] = driftp;
  return _dist->r(params, rng);
}

} //namespace wiender
} //namespace jags
