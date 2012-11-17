#include <config.h>
#include "DWiener.h"

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

namespace wiener {

DWiener::DWiener() : ScalarDist("dwiener", 4, DIST_UNBOUNDED)
{}

double DWiener::logDensity(double x, PDFType type,
			vector<double const *> const &parameters,
			double const *lower, double const *upper) const
{
  if (lower && x < *lower) return JAGS_NEGINF;
  if (upper && x > *upper) return JAGS_NEGINF;
  if (upper && lower && *upper < *lower) return JAGS_NEGINF;
    
  double loglik =  d(x, type, parameters, true);

  if (lower || upper) {
    /* In theory, we just have to subtract log[P(lower <= X <=
             upper)] from the log likelihood. But we need to work around
             numerical problems. */
    bool have_lower = lower && p(*lower, parameters, true, false) > 0;
    bool have_upper = upper && p(*upper, parameters, false, false) > 0;

    if (have_lower && have_upper) {
      if (p(*lower, parameters, false, false) < 0.5) {
        //Use upper tail
        loglik -= log(p(*lower, parameters, false, false) -
          p(*upper, parameters, false, false));
      }
      else {
        //Use lower tail
        loglik -= log(p(*upper, parameters, true, false) - 
          p(*lower, parameters, true, false));
      }
    }
    else if (have_lower) {
      loglik -= p(*lower, parameters, false, true);
    }
    else if (have_upper) {
      loglik -= p(*upper, parameters, true, true);
    }
  }
  return loglik;
}


double DWiener::randomSample(vector<double const *> const &parameters,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
  /*
  if (lower || upper) {
    double plower = 0, pupper = 1;
    if (lower) {
        plower = p(*lower, parameters, true, false);
    }
    if (upper) {
        pupper = p(*upper, parameters, true, false);
    }
    double u = plower + rng->uniform() * (pupper - plower);
    return q(u, parameters, true, false);
  }
  else {
    return r(parameters, rng);
  }
  */
  return r(parameters, rng);

}

double DWiener::typicalValue(vector<double const *> const &parameters,
			  double const *lower, double const *upper) const
{
  double llimit = l(parameters), ulimit = u(parameters);
  double plower = 0, pupper = 1;
  
  if (lower) {
    llimit = max(llimit, *lower);
    plower = p(llimit, parameters, true, false);
  }

  if (upper) {
    ulimit = min(ulimit, *upper);
    pupper = p(ulimit, parameters, true, false);
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
            BIAS(par)<1 && BIAS(par)>0 &&
            TER(par) > 0);
}

double DWiener::d(double x, PDFType type, 
         vector<double const *> const &par, bool give_log) const
{
  double v=DRIFT(par), w=BIAS(par);
  double kl, ks, ans;
  int k,K;
  //int acc; // accuracy - true if response was correct 

#ifdef IEEE_754
  if (jags_isnan(x) || jags_isnan(v) || jags_isnan(BOUND(par)) || jags_isnan(w) || jags_isnan(TER(par)))
    return x + v + BOUND(par) + w + TER(par);
#endif
  if (!jags_finite(x) || !jags_finite(BOUND(par))) return (give_log ? JAGS_NEGINF : 0);
  if (w < 0 || w > 1 || BOUND(par) <= 0 || TER(par) <= 0) return JAGS_NAN;

  // extract RT and accuracy from x
  if (x<0) {
    //acc = false;
    x = fabs(x);
  }
  else {
    //acc = true;
    w = 1-w;
    v = -v;
  }
  
  x = x-TER(par); // remove non-decision time from x
  x = x/pow(BOUND(par),2); // convert t to normalized time tt

  // calculate number of terms needed for large t
  if (M_PI*x*WIENER_ERR<1) { // if error threshold is set low enough
      kl=sqrt(-2*log(M_PI*x*WIENER_ERR)/(pow(M_PI,2)*x)); // bound
      kl=(kl>1/(M_PI*sqrt(x))) ? kl : 1/(M_PI*sqrt(x)); // ensure boundary conditions met
  }
  else { // if error threshold set too high
      kl=1/(M_PI*sqrt(x)); // set to boundary condition
  }
  // calculate number of terms needed for small t
  if ((2*sqrt(2*M_PI*x)*WIENER_ERR)<1) { // if error threshold is set low enough
      ks=2+sqrt(-2*x*log(2*sqrt(2*M_PI*x)*WIENER_ERR)); // bound
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
    ans+((-v*BOUND(par)*w -(pow(v,2))*(x*pow(BOUND(par),2))/2)-log(pow(BOUND(par),2))) : 
    ans*exp(-v*BOUND(par)*w -(pow(v,2))*(x*pow(BOUND(par),2))/2)/(pow(BOUND(par),2));
}

double DWiener::p(double q, vector<double const *> const &par, bool lower, 
    bool give_log) const
{
  // TODO
  return JAGS_NAN;
}

double DWiener::q(double p, vector<double const *> const &par, bool lower, 
    bool log_p) const
{
  // TODO
  if(p == 0.5) return 1;
  else return JAGS_NAN;
}

double DWiener::r(vector<double const *> const &par, RNG *rng) const
{
  r_random_walk(par,rng); // if dt is to small, it won't work!
}

double DWiener::r_random_walk(vector<double const *> const &par, RNG *rng, double dt) const
{
  double t,sigma=1;
	unsigned int i = 0;
	double y = BIAS(par)*BOUND(par);
	while(y < BOUND(par) && y > 0) {
		y = y + sqrt(dt) * (sigma*rng->normal()+DRIFT(par));
		i+=1;
	}
  if(y >= BOUND(par)) t = -(i*dt+TER(par));
  else t = i*dt+TER(par);
	return t;
}

double DWiener::r_rejection_based(vector<double const *> const &par, RNG *rng) const
{
  double t=0;
  vector<double > tmp_pars(2,0);
  tmp_pars[0] = *par[0];
  tmp_pars[1] = *par[3];
  if(BIAS(par) == 0.5) {
    t = r_rejection_based_symmetric(tmp_pars,rng);
  }
  else {
    double tmp_t,a,z= BIAS(par)*BOUND(par);
    bool bound_hit=false;
    while (!bound_hit) {
      // near upper bound
      if (z/BOUND(par) > .5) {
        a = (BOUND(par)-z)*2;
        tmp_pars[0] = a;
        tmp_t = r_rejection_based_symmetric(tmp_pars,rng);
        t += std::abs(tmp_t);
        if(tmp_t < 0) {
          bound_hit = true;
          t = -t;
        }
        else z=BOUND(par)-a;
      }
      // near lower bound
      else if (z/BOUND(par) < .5) {
        a = z*2;
        tmp_pars[0] = a;
        tmp_t = r_rejection_based_symmetric(tmp_pars,rng);
        t += std::abs(tmp_t);
        if(tmp_t < 0) z = z*2;
        else bound_hit=true;
      }
      // symmetric (in the middle of both bounds)
      else {
        a = (BOUND(par)-z)*2;
        tmp_pars[0] = a;
        tmp_t = r_rejection_based_symmetric(tmp_pars,rng);
        t += std::abs(tmp_t);
        if(tmp_t < 0) {
          bound_hit = true;
          t = -t;
        }
        else bound_hit = true;
      }
    } // end while
  }
  if (t >= 0)  t = (t+TER(par));
  else t = (t-TER(par));
  return t;
}

double DWiener::r_rejection_based_symmetric(vector<double> par, RNG *rng) const
{
  // Rejection based algorithm - see Tuerlinckx et al. (2001) 
  double u,t,crit;
  bool converged=false;
  double mu = par[1]; 
  double sigma = 1;
  double a = par[0];
  double lambda = pow(mu,2)/(2*pow(sigma,2)) + (pow(M_PI,2)*pow(sigma,2))/(2*pow(a,2));
  double M = ( (M_PI*pow(sigma,2))/pow(a,2) 
      * ( exp((a*mu)/(2*pow(sigma,2))) + exp((-a*mu)/(2*pow(sigma,2))) ) ) 
      / (pow(mu,2)/(2*pow(sigma,2)) + (pow(M_PI,2)*pow(sigma,2))/(2*pow(a,2)));
  double F = (pow(M_PI,2)*pow(sigma,4)) / (pow(M_PI,2)*pow(sigma,4)+pow(mu,2)*pow(a,2));
  double infinite_sum, infinite_sum_o1=0, infinite_sum_o2=0;

  do {
    u = rng->uniform();
    t = -1/lambda * std::abs(log(1-u));
    converged = false;
    infinite_sum=0;
    for (int n=1;!converged;n++) {
      infinite_sum_o2 = infinite_sum_o1;
      infinite_sum_o1 = infinite_sum;
      infinite_sum += (2*n+1)*pow(-1,n)*pow(1-u,F*pow(2*n+1,2));
      /* stop approximation of ininite_sum,
       * when the difference to the (n-1) sum and the (n-2) sum is below
       * WIENER_ERR. As proposed by Tuerlinckx (2004). */
      if (n>2 && std::abs(infinite_sum-infinite_sum_o1) < WIENER_ERR && std::abs(infinite_sum-infinite_sum_o2) < WIENER_ERR) converged=true;
      else if (n>10000) break; // break for-loop if it does not converge
    }
    if(converged) crit = 1 + pow((1-u),F) * infinite_sum;
    else crit = 0; // if not converged, try a different value.
  } while((rng->uniform())>crit); // take the t, if within target distribution

  if (((rng->normal()+mu)*sigma)>=0) return -t; // hit the upper bound
  else return t; // hit the lower bound
}

}
