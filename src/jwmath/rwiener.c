/* 
 *  rwiener function
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
 *    double rwiener(double v, double a, double w, double ter, RNG *rng)
 *
 *  DESCRIPTION
 *
 *    Random variates from the diffusion model distribuion.
 *    Based on code by Francis Tuerlinckx.
 *
 */

#include "jwmath.h"

/* TODO: currently put on ice, other things are more important...
double rwiener_new(double v, double a, double w, double ter, RNG *rng)
{
  double sigma=1,sample=0,allsample=0,lambda,radius;
  int repeat=1;

  if (w == 0.5) {
    lambda = ( pow(v,2)/(2*pow(sigma,2)) + (pow(M_PI,2)*pow(sigma,2))/(2*pow(a,2)) )
    prob = 0.5 // todo replace with proper value
    if (unif_rand(rng) > prob) { // X=1, right answear, upper boundary
      do {
        sample = lambda * exp_rand(rng);
        if ( unif_rand(rng)<= dwiener(-sample,v,a,w,ter,0) / (exp(-sample / lambda) / lambda) ) {
          repeat = 0;
        }
      } while(repeat);
    }
    else { // X=0, wrong answear, lower boundary
      do {
        sample = lambda * exp_rand(rng);
        if ( unif_rand(rng)<= dwiener(sample,v,a,w,ter,0) / (exp(-sample / lambda) / lambda) ) {
          repeat = 0;
        }
      } while(repeat);
    }
  }
  else { // w != 0.5
    lambda = 0;
    prop = 0.5;
    radius = ((a-w*a)>(w*a)?(a-w*a):(w*a));
    do {
      if (unif_rand(rng) > prob) { // X=1, right answear, upper boundary
        do {
          sample = lambda * exp_rand(rng);
          if ( unif_rand(rng)<= dwiener(-sample,v,radius,w,ter,0) / (exp(-sample / lambda) / lambda) ) {
            repeat = 0;
            radius = radius + radius/2;
          }
        } while(repeat);
      }
      else { // X=0, wrong answear, lower boundary
        do {
          sample = lambda * exp_rand(rng);
          if ( unif_rand(rng)<= dwiener(sample,v,radius,w,ter,0) / (exp(-sample / lambda) / lambda) ) {
            repeat = 0;
            radius = radius - radius/2;
          }
        } while(repeat);
        allsample += sample
        if (radius >= a) sample -= allsample;
        else if (radius <= 0) sample = allsample;
        else repeat =1;
      }
    } while(repeat);
  }
  return sample;
} 
*/ 

double rwiener(double v, double a, double w, double ter, RNG *rng)
{
  double D=.005,totaltime,startpos,radius,lambda,F,prob,tt,dir_,l,s1,s2,tnew,t_delta;
  int uu,i=0;
  
  a/=10;
  w*=a; // convert w to z
  v/=10;
      
  totaltime=0;
  startpos=0;
  radius=(fabs(a-w)>fabs(-w))?fabs(a-w):fabs(-w);
  
  while(1)
  {
      if (v==0) {
          lambda = 0.25*D*pow(M_PI,2)/pow(radius,2);
          F=1;
          prob = .5;
      } 
      else {
          lambda = 0.25*v*v/D + 0.25*D*pow(M_PI,2)/pow(radius,2);
          F=D*M_PI/(radius*v);
          F=F*F/(1+F*F);
          prob=exp(radius*v/D);
          prob=prob/(1+prob);
      }
      dir_= unif_rand(rng)<prob ? 1 : -1;
      
      do {
          s2 = unif_rand(rng);
          s1 = unif_rand(rng);
          tnew=0;
          t_delta=0;
          uu=0;
          while( (fabs(t_delta)>0) | (!uu) )
          {
              tt = 2*++uu+1;
              t_delta = tt * (uu%2?-1:1) * pow(s1,(F*pow(tt,2)));
              tnew += t_delta;
          }
          
          l = 1 + pow(s1,-F) * tnew;
      } while (s2>l);
      
      totaltime+=fabs(log(s1))/lambda;
      dir_=startpos+dir_*radius;
      
      if (dir_>=(a-w))
      {
          return (totaltime+ter); // acc=1, therefore +RT;
      }
      else
      {
          if (dir_<=(-w))
          {
              return -(totaltime+ter); // acc=0, therefore -RT;
          }
          else
          {
              startpos=dir_;
              radius=(fabs((a-w)-startpos)>fabs((-w)-startpos))?fabs((a-w)-startpos):fabs((-w)-startpos);
          }
      }
      i++;
      if(i>499) return JWML_NAN;
  }
}
