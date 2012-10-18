#include <config.h>
#include "DWiener.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>

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
  double v=DRIFT(par), a=BOUND(par), w=RELST(par), ter=TER(par);
  double kl, ks, ans;
  int k,K;
  int acc; // accuracy - true if response was correct
  // TODO: Do we need this variable anyway?

#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(v) || ISNAN(a) || ISNAN(w) || ISNAN(ter))
    return x + v + a + w + ter;
#endif
  if (!R_FINITE(x) || !R_FINITE(a)) return (give_log ? JWML_NEGINF : 0);
  if (w < 0 || w > 1 || a <= 0 || ter <= 0) return JWML_NAN;

  // extract RT and accuracy from x
  if (x<0) {
    acc = 0; // false
    x = fabs(x);
  }
  else {
    acc = 1; // true
    w = 1-w;
    v = -v;
  }
  
  x = x-ter; // remove non-decision time from x
  x = x/pow(a,2); // convert t to normalized time tt

  // calculate number of terms needed for large t
  if (M_PI*x*JWML_WIENER_ERR<1) { // if error threshold is set low enough
      kl=sqrt(-2*log(M_PI*x*JWML_WIENER_ERR)/(pow(M_PI,2)*x)); // bound
      kl=(kl>1/(M_PI*sqrt(x))) ? kl : 1/(M_PI*sqrt(x)); // ensure boundary conditions met
  }
  else { // if error threshold set too high
      kl=1/(M_PI*sqrt(x)); // set to boundary condition
  }
  // calculate number of terms needed for small t
  if ((2*sqrt(2*M_PI*x)*JWML_WIENER_ERR)<1) { // if error threshold is set low enough
      ks=2+sqrt(-2*x*log(2*sqrt(2*M_PI*x)*JWML_WIENER_ERR)); // bound
      ks=(ks>sqrt(x)+1) ? ks : sqrt(x)+1; // ensure boundary conditions are met
  }
  else { // if error threshold was set too high
      ks=2; // minimal kappa for that case
  }

  // compute density: f(tt|0,1,w)
  ans=0; //initialize density
  if (ks<kl) { // if small t is better (i.e., lambda<0)
      K=ceil(ks); // round to smallest integer meeting error
      for (k=-floor((K-1)/2); k<=ceil((K-1)/2); k++) { // loop over k
          ans=ans+(w+2*k)*exp(-(pow((w+2*k),2))/2/x); // increment sum
      }
      ans= give_log ? log(ans)-0.5*log(2)-M_LN_SQRT_PI-1.5*log(x) : ans/sqrt(2*M_PI*pow(x,3)); // add constant term
  }
  else { // if large t is better...
      K=ceil(kl); // round to smallest integer meeting error
      for (k=1; k<=K; k++) {
          ans=ans+k*exp(-(pow(k,2))*(pow(M_PI,2))*x/2)*sin(k*M_PI*w); // increment sum
      }
      ans= give_log ? log(ans)+2*M_LN_SQRT_PI : ans*M_PI; // add constant term
  }

  // convert to f(t|v,a,w) and return result
  return give_log ? 
    ans+((-v*a*w -(pow(v,2))*(x*pow(a,2))/2)-log(pow(a,2))) : 
    ans*exp(-v*a*w -(pow(v,2))*(x*pow(a,2))/2)/(pow(a,2));
}

double
DWiener::p(double q, vector<double const *> const &par, bool lower, bool give_log)
  const
{
    return JWML_NAN;
}

double 
DWiener::q(double p, vector<double const *> const &par, bool lower, bool log_p)
  const
{
    return JWML_NAN;
}

double 
DWiener::r(vector<double const *> const &par, RNG *rng) const
{
    return JWML_NAN;
}

}
