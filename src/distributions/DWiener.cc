#include <config.h>
#include "DWiener.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <cmath>
#include <algorithm>

#include <jwmath/jwmath.h>

using std::string;
using std::vector;
using std::log;
using std::min;
using std::max;

#define BOUND(par) (*par[0])
#define TER(par) (*par[1])
#define RELST(par) (*par[2])
#define DRIFT(par) (*par[3])

namespace wiener {

DWiener::DWiener()
    : ScalarDist("dwiener", 4, DIST_UNBOUNDED)
{}

double 
DWiener::logDensity(double x, PDFType type,
			vector<double const *> const &parameters,
			double const *lower, double const *upper) const
{
    if (lower && x < *lower)
	return JAGS_NEGINF;
    if (upper && x > *upper)
	return JAGS_NEGINF;
    if (upper && lower && *upper < *lower)
	return JAGS_NEGINF;
    
    double loglik =  d(x, type, parameters, true);

    if (type != PDF_PRIOR && (lower || upper)) {
	//Normalize truncated distributions

	double ll = 0;
	if (lower) {
	    ll = *lower;
	}

	/* In theory, we just have to subtract log[P(lower <= X <=
           upper)] from the log likelihood. But we need to work around
           numerical problems. */

	bool have_lower = lower && p(ll, parameters, true, false) > 0;
	bool have_upper = upper && p(*upper, parameters, false, false) > 0;

	if (have_lower && have_upper) {
	    if (p(ll, parameters, false, false) < 0.5) {
		//Use upper tail
		loglik -= log(p(ll, parameters, false, false) -
			      p(*upper, parameters, false, false));
	    }
	    else {
		//Use lower tail
		loglik -= log(p(*upper, parameters, true, false) - 
			      p(ll, parameters, true, false));
	    }
	}
	else if (have_lower) {
	    loglik -= p(ll, parameters, false, true);
	}
	else if (have_upper) {
	    loglik -= p(*upper, parameters, true, true);
	}
    }

    return loglik;
}

double 
DWiener::randomSample(vector<double const *> const &parameters,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    if (lower || upper) {

	double plower = 0, pupper = 1;
	if (lower) {
	    plower = calPlower(*lower, parameters);
	}
	if (upper) {
	    pupper = calPupper(*upper, parameters);
	}
	
	double u = plower + rng->uniform() * (pupper - plower);
	return q(u, parameters, true, false);
    }
    else {
	return r(parameters, rng);
    }

}

double DWiener::calPlower(double lower, 
			      vector<double const*> const &parameters) const
{
  //P(X < lower)
	return p(lower, parameters, true, false);
}

double DWiener::calPupper(double upper,
			      vector<double const*> const &parameters) const
{
  //P(X <= upper)
  return p(upper, parameters, true, false);
}

double 
DWiener::typicalValue(vector<double const *> const &parameters,
			  double const *lower, double const *upper) const
{
    double llimit = l(parameters), ulimit = u(parameters);
    double plower = 0, pupper = 1;
    
    if (lower) {
	llimit = max(llimit, *lower);
	plower = calPlower(llimit, parameters);
    }

    if (upper) {
	ulimit = min(ulimit, *upper);
	pupper = calPupper(ulimit, parameters);
    }
    
    double pmed = (plower + pupper)/2;
    double med = q(pmed, parameters, true, false);	

    //Calculate the log densities
    double dllimit = d(llimit, PDF_FULL, parameters, true);
    double dulimit = d(ulimit, PDF_FULL, parameters, true);
    double dmed = d(med, PDF_FULL, parameters, true);

    //Pick the median if it has the highest density, otherwise pick
    //a point near to (but not on) the boundary
    if (dmed >= dllimit && dmed >= dulimit) {
	return med;
    }
    else if (dulimit > dllimit) {
	return q(0.1 * plower + 0.9 * pupper, parameters, true, false);
    }
    else {
	return q(0.9 * plower + 0.1 * pupper, parameters, true, false);
    }
}

bool DWiener::checkParameterValue (vector<double const *> const &par) const
{
    return (BOUND(par)>0 && 
        (RELST(par)<1 && RELST(par)>0) &&
        TER(par) > 0);
}

double
DWiener::d(double x, PDFType type, 
         vector<double const *> const &par, bool give_log) const
{
    return dwiener(x, DRIFT(par), BOUND(par), RELST(par), TER(par), give_log);
}

double
DWiener::p(double q, vector<double const *> const &par, bool lower, bool give_log)
  const
{
    return pwiener(q, DRIFT(par), BOUND(par), RELST(par), TER(par), lower, give_log);
}

double 
DWiener::q(double p, vector<double const *> const &par, bool lower, bool log_p)
  const
{
    return qwiener(p, DRIFT(par), BOUND(par), RELST(par), TER(par), lower, log_p);
}

double 
DWiener::r(vector<double const *> const &par, RNG *rng) const
{
    return rwiener(DRIFT(par), BOUND(par), RELST(par), TER(par), rng);
}

}
