/* 
 *  dwiener function
 *  Copyright (C) 2012 Dominik Wabersich <dominik.wabersich@gmail.com>
 *  and Joachim Vandekerckhove <joachim@uci.edu>
 *
 *  When using these functions, please cite as: 
 *      Wabersich, D. and Vandekerckhove, J. (in preparation). Hacking JAGS: 
 *      Adding custom distributions to JAGS, with a diffusion model example.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 *  SYNOPSIS
 *    
 *    #include "JRmath.h"
 *    double dwiener(double x, double v, double a, double w, double ter, int give_log)
 *
 *  DESCRIPTION
 *
 *    Calculates the density for a given value of x and given the
 *    parameters. 
 *    Based on calculations by Daniel Navarro & Ian Fuss.
 *    x contains both the accuracy and RT variable,
 *    by storing the binary accuracy through positive or negative RT value.
 *    RT is therefore abs(RT).
 *
 */

#include "jwmath.h"

double dwiener(double x, double v, double a, double w, double ter, int give_log)
{
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
