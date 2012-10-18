/* 
 *  qwiener function
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
 *    double qwiener(double p, double v, double a, double w, double ter, int lower, int log_p)
 *
 *  DESCRIPTION
 *
 *    Calculates the quantile for a given value x.
 *
 */

#include "jwmath.h"

double qwiener(double p, double v, double a, double w, double ter, int lower, int log_p)
{
  return JWML_NAN;
}

