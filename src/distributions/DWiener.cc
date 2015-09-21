#include <config.h>
#include "DWiener.h"

#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>
#include <JRmath.h>

using std::vector;
using std::log;
using std::min;
using std::max;
using std::string;

#define BOUND(par) (*par[0])
#define TER(par) (*par[1])
#define BIAS(par) (*par[2])
#define DRIFT(par) (*par[3])

namespace jags {
namespace wiener {

DWiener::DWiener() : ScalarDist("dwiener", 4, DIST_UNBOUNDED)
{}

DWiener::DWiener(string const &name, unsigned int npar) : ScalarDist(name, npar, DIST_UNBOUNDED)
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
    bool have_lower = lower && p_full(*lower, parameters, true, false) > 0;
    bool have_upper = upper && p_full(*upper, parameters, false, false) > 0;

    if (have_lower && have_upper) {
      if (p(*lower, parameters, false, false) < 0.5) {
        //Use upper tail
        loglik -= log(p(*lower, parameters, false, false) -
          p_full(*upper, parameters, false, false));
      }
      else {
        //Use lower tail
        loglik -= log(p(*upper, parameters, true, false) - 
          p_full(*lower, parameters, true, false));
      }
    }
    else if (have_lower) {
      loglik -= p_full(*lower, parameters, false, true);
    }
    else if (have_upper) {
      loglik -= p_full(*upper, parameters, true, true);
    }
  }
  return loglik;
}


double DWiener::randomSample(vector<double const *> const &parameters,
        double const *lower, double const *upper, RNG *rng) const
{
  if (lower || upper) {
    double plower = 0, pupper = 1;
    if (lower) {
        plower = p_full(*lower, parameters, true, false);
    }
    if (upper) {
        pupper = p_full(*upper, parameters, true, false);
    }
    double u = plower + rng->uniform() * (pupper - plower);
    return q_full(u, parameters, true, false);
  }
  else {
    return r(parameters, rng);
  }

}

double DWiener::typicalValue(vector<double const *> const &parameters,
        double const *lower, double const *upper) const
{
  double llimit = l(parameters), ulimit = u(parameters);
  double plower = 0, pupper = 1;
  
  if (lower) {
    llimit = max(llimit, *lower);
    plower = p_full(llimit, parameters, true, false);
  }

  if (upper) {
    ulimit = min(ulimit, *upper);
    pupper = p_full(ulimit, parameters, true, false);
  }
  
  double pmed = (plower + pupper)/2;
  double med = q_full(pmed, parameters, true, false);  

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
    return q_full(0.1 * plower + 0.9 * pupper, parameters, true, false);
  }
  else {
    return q_full(0.9 * plower + 0.1 * pupper, parameters, true, false);
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

/* This function takes a quantil and returns the full probability
 * of hitting any of the two boundaries */
double DWiener::p_full(double q, vector<double const *> const &par, bool lower, 
    bool give_log) const
{
  double ptmp;
  if (q < 0) return JAGS_NAN;
  if(q == JAGS_POSINF) return JAGS_POSINF;
  ptmp = (p(q,par,true,false) + p(-q,par,true,false));
  if (!lower) return give_log?log(1-ptmp):(1-ptmp);
  else return give_log?log(ptmp):ptmp;
}

double DWiener::p(double q, vector<double const *> const &par, bool lower, 
    bool give_log) const
{
  // TODO: if not lower
  if(!lower) return JAGS_NAN;
  if(q == JAGS_POSINF || q == JAGS_NEGINF) return JAGS_POSINF;
  double p;
  if (jags_isnan(q)) return JAGS_NAN;
  if (fabs(q) <= TER(par)) return give_log?JAGS_NEGINF:0;
  if (q < 0) { // lower boundary 0
    p = F_lower(fabs(q)-TER(par), DRIFT(par), BOUND(par), BIAS(par));
  }
  else { // upper boundary a
    p = F_lower(q-TER(par), (-DRIFT(par)), BOUND(par), (1-BIAS(par)));
  }
  // TODO: Make calculations more efficient by using give_log
  if (give_log) return log(p);
  else return p;
}

double DWiener::F_lower(double q, double v, double a, double w) const
{
  /*  double sigma = 1;
      a = a / sigma;
      v = v / sigma; */
  double F;
  int K_l = K_large(q, v, a, w);
  int K_s = K_small(q, v, a, w);
  if (K_l < 10*K_s) F = Fl_lower(q, v, a, w, K_l);
  else F = Fs_lower(q, v, a, w, K_s);
  return F;
}

double DWiener::Fl_lower(double q, double v, double a, double w, int K) const
{
  double F=0;
  for(int k=K; k>=1; k--) F = F - k / (v*v*1.0 + k*k*M_PI*M_PI/(a*1.0)/a) * exp(-v*a*w*1.0 - 0.5*v*v*q - 0.5*k*k*M_PI*M_PI/(a*1.0)/a*q) * sin(M_PI*k*w);
  return prob_upperbound(v, a, w) + 2.0*M_PI/(a*1.0)/a * F;
}

double DWiener::Fs_lower(double q, double v, double a, double w, int K) const
{
  if (v == 0) return(Fs0_lower(q, a, w, K));
  double S1=0,S2=0;
  double sqt = sqrt(q);
  for(int k=K; k>=1; k--) {
    S1 = S1 + exp_pnorm(2*v*a*k, -sign(v)*(2*a*k+a*w+v*q)/sqt) -
           exp_pnorm(-2*v*a*k - 2*v*a*w, sign(v)*(2*a*k+a*w-v*q)/sqt);
    S2 = S2 + exp_pnorm(-2*v*a*k, sign(v)*(2*a*k-a*w-v*q)/sqt) -
           exp_pnorm(2*v*a*k - 2*v*a*w, -sign(v)*(2*a*k-a*w+v*q)/sqt);
  }
  return prob_upperbound(v, a, w) + sign(v) * ((pnorm(-sign(v) * (a*w+v*q)/sqt,0,1,1,0) -
           exp_pnorm(-2*v*a*w, sign(v) * (a*w-v*q)/sqt)) + S1 + S2);
}

double DWiener::Fs0_lower(double q, double a, double w, int K) const
{
  double F=0;
  for(int k=K; k>=0; k--) {
    F = F - pnorm((-2*k - 2 + w)*a/sqrt(q),0,1,1,0) + pnorm((-2*k - w)*a/sqrt(q),0,1,1,0);
  }

  return 2*F;
}

double DWiener::prob_upperbound(double v, double a, double w) const
{
  double e = exp(-2.0 * v * a * (1.0-w));
  if(e == JAGS_POSINF) return 1;
  else if(v == 0 || w == 1) return (1-w);
  else return ((1 - e) / (exp(2.0*v*a*w) - e));
}

double DWiener::exp_pnorm(double a, double b) const
{
  double r;
  if (jags_isnan(r) && b < -5.5) r = 1/sqrt(2) * exp(a - b*b/2) * (0.5641882/b/b/b - 1/b/sqrt(M_PI));
  else r = exp(a) * pnorm(b,0,1,1,0);
  return r;
}

int DWiener::K_large(double q, double v, double a, double w)  const
{
  double sqrtL1 = sqrt(1/q) * a/M_PI;
  double sqrtL2 = sqrt(max(1.0, -2/q*a*a/M_PI/M_PI * (log(WIENER_ERR*M_PI*q/2 * (v*v + M_PI*M_PI/a/a)) + v*a*w + v*v*q/2)));
  return ceil(max(sqrtL1, sqrtL2));
}

int DWiener::K_small(double q, double v, double a, double w, double epsilon)  const
{
  if(v == 0) return ceil(max(0.0, w/2 - sqrt(q)/2/a * qnorm(max(0.0, min(1.0, epsilon/(2-2*w))),0,1,1,0)));
  if(v > 0) return(K_small(q, -v, a, w, exp(-2*a*w*v)*epsilon));
  double S2 = w - 1 + 0.5/v/a * log(epsilon/2 * (1-exp(2*v*a)));
  double S3 = (0.535 * sqrt(2*q) + v*q + a*w)/2/a;
  double S4 = w/2 - sqrt(q)/2/a * qnorm(max(0.0, min(1.0, epsilon * a / 0.3 / sqrt(2*M_PI*q) * exp(v*v*q/2 + v*a*w))),0,1,1,0);
  return ceil(max(max(max(S2, S3), S4), 0.0));
}

/* This function takes the full probability of hitting any boundary
 * and returns the appropriate quantile */
double DWiener::q_full(double p, vector<double const *> const &par, bool lower, 
    bool log_p) const
{
  double p_tmp,pmin,pmax,pmid;
  double qmin,qmax,q;

  if (log_p) p_tmp = exp(p);
  else p_tmp = p;
  if (p_tmp > 1) return JAGS_NAN;
  if(!lower) p_tmp = 1-p_tmp;

  q = 1;
  qmin = 0;
  pmin = 0;
  qmax = JAGS_POSINF;
  pmax = 1;

  int c=0;
  do {
    c++;
    pmid = this->p_full(q, par, true, false);
    if (fabs(p_tmp)<=pmid) { // near lower point
      pmax = pmid;
      qmax = q;
      q = qmin + (qmax-qmin)/2;
    }
    else { // near upper point
      pmin = pmid;
      qmin = q;
      if (!(qmax == JAGS_NEGINF || qmax == JAGS_POSINF)) 
        q = qmin + (qmax-qmin)/2;
      else
        q = q*10;
    }
    if(jags_isnan(pmid)) return JAGS_NAN;
    if(q>=1e+10) return JAGS_POSINF;
  } while(fabs(p_tmp-pmid) > WIENER_ERR && c < 1000); // defines the accuracy

  return q;
}
double DWiener::q(double p, vector<double const *> const &par, bool lower, 
    bool log_p) const
{
  double p_tmp,pmin,pmax,pmid;
  double qmin,qmax,q;

  if (log_p) p_tmp = exp(p);
  else p_tmp = p;
  if (fabs(p_tmp) > 1) return JAGS_NAN;

  q = 1;
  qmin = 0;
  pmin = 0;
  qmax = JAGS_POSINF;
  pmax = 1;

  int c=0;
  do {
    c++;
    if (p_tmp>=0) pmid = this->p(q, par, lower, false);
    else pmid = this->p(-q, par, lower, false);
    if (fabs(p_tmp)<=pmid) { // near lower point
      pmax = pmid;
      qmax = q;
      q = qmin + (qmax-qmin)/2;
    }
    else { // near upper point
      pmin = pmid;
      qmin = q;
      if (!(qmax == JAGS_NEGINF || qmax == JAGS_POSINF)) 
        q = qmin + (qmax-qmin)/2;
      else
        q = q*10;
    }
    if(jags_isnan(pmid)) return JAGS_NAN;
    if(q>=1e+10) return JAGS_POSINF;
  } while(fabs(fabs(p_tmp)-pmid) > WIENER_ERR && c < 1000); // defines the accuracy

  return q;
}

double DWiener::r(vector<double const *> const &par, RNG *rng) const
{
  //return r_rejection_based(par,rng);
  //return r_random_walk(par,rng); // if dt is to small, it won't work!
  return r_rejection_based_2(BOUND(par), TER(par), BOUND(par)*BIAS(par), DRIFT(par),rng);
}

/*
double DWiener::r_random_walk(vector<double const *> const &par, RNG *rng, double dt) const
{
  double t,sigma=1;
  double p = .5 * (1+((DRIFT(par)*sqrt(dt))/sigma));
  //double q = .5 * (1-((mu*sqrt(dt))/sigma));
  int i = 0;
  double y = BIAS(par)*BOUND(par);
  while(y < BOUND(par) && y > 0)
  {
    if(rng->uniform() <= p) y = y + sigma*sqrt(dt);
    else y = y - sigma*sqrt(dt);
    i++;
  }
  if(y >= BOUND(par)) t = (i*dt+TER(par));
  else t = -(i*dt+TER(par));
  return t;

}

// NOTE: this random walk does not work!
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

// TODO: Need to find error in this function...
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
        t += fabs(tmp_t);
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
        t += fabs(tmp_t);
        if(tmp_t < 0) z = z*2;
        else bound_hit=true;
      }
      // symmetric (in the middle of both bounds)
      else {
        a = (BOUND(par)-z)*2;
        tmp_pars[0] = a;
        tmp_t = r_rejection_based_symmetric(tmp_pars,rng);
        t += fabs(tmp_t);
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
    t = -1/lambda * fabs(log(1-u));
    converged = false;
    infinite_sum=0;
    for (int n=1;!converged;n++) {
      infinite_sum_o2 = infinite_sum_o1;
      infinite_sum_o1 = infinite_sum;
      infinite_sum += (2*n+1)*pow(-1,n)*pow(1-u,F*pow(2*n+1,2));
      // stop approximation of ininite_sum,
      // when the difference to the (n-1) sum and the (n-2) sum is below
      // WIENER_ERR. As proposed by Tuerlinckx (2004).
      if (n>2 && fabs(infinite_sum-infinite_sum_o1) < WIENER_ERR && fabs(infinite_sum-infinite_sum_o2) < WIENER_ERR) converged=true;
      else if (n>10000) break; // break for-loop if it does not converge
    }
    if(converged) crit = 1 + pow((1-u),F) * infinite_sum;
    else crit = 0; // if not converged, try a different value.
  } while((rng->uniform())>crit); // take the t, if within target distribution

  if (((rng->normal()+mu)*sigma)>=0) return -t; // hit the upper bound
  else return t; // hit the lower bound
}
*/

double DWiener::r_rejection_based_2(double a, double ter, double z, double v, RNG *rng) const
{
  /*  mere copy of wdmrnd.cpp by JV, only changes:
   *  - return value double instead of void
   *  - removed *t and *x, instead returning t or -t 
   *  - added variable t (double)
   *  - replaced GNU gsl with JAGS rng
   *  - absol replaced with fabs
   *  - amin replaced with min
   *  - pi replaced with M_PI
   */
  double dt=1e-15,tau=.1,D=.005,totaltime,startpos,ndrt,
  zz,Aupper,Alower,radius,lambda,F,prob,tt,dir_,l,s1,s2,tnew,t_delta;
  int uu,i;
  bool finish;
  double t;
  
  a/=10;
  z/=10;
  v/=10;

  finish = false;
  totaltime=0;
  startpos=0;
  Aupper=a-z;
  Alower=-z;
  radius=min(fabs(Aupper),fabs(Alower));
  
  while (!finish) {
    if (v==0){
      lambda = 0.25*D*M_PI*M_PI/(radius*radius);
      F=1;
      prob = .5;
    } 
    else {
      lambda = 0.25*v*v/D + 0.25*D*M_PI*M_PI/(radius*radius);
      F=D*M_PI/(radius*v);
      F=F*F/(1+F*F);
      prob=exp(radius*v/D);
      prob=prob/(1+prob);
    }
    dir_= rng->uniform()<prob ? 1 : -1;
    l=-1;
    s2=0;
    
    while (s2>l) {
      s2 = rng->uniform();
      s1 = rng->uniform();
      tnew=0;
      t_delta=0;
      uu=0;
      
      while ( (fabs(t_delta)>dt) | (!uu) ) {
        tt = 2*++uu+1;
        t_delta = tt * (uu%2?-1:1) * pow(s1,(F*tt*tt));
        tnew += t_delta;
      }
      
      l = 1 + pow(s1,-F) * tnew;
    }/*end while (s2>l) */
    
    totaltime+=fabs(log(s1))/lambda;
    dir_=startpos+dir_*radius;
    
    if (dir_+dt>Aupper) {
      //*t=totaltime+ter;
      //*x=1;
      t = totaltime+ter;
      finish=true;
      return t;
    }
    else {
      if (dir_-dt<Alower) {
        //*t=totaltime+ter;
        //*x=0;
        t = -(totaltime+ter);
        finish=true;
        return t; 
      }
      else {
        startpos=dir_;
        radius=min(fabs(Aupper-startpos),fabs(Alower-startpos));
      }
    }
  } /*end while (!finish) */
}

} //namespace wiener
} //namespace jags
