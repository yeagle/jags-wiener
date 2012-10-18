/* JAGS Wiener Mathlib */

/*
 *  Copyright (C) 2012 Dominik Wabersich <dominik.wabersich@gmail.com>
 *  and Joachim Vandekerckhove <joachim@uci.edu>
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
 */

// to resolve unnecessary R lib dependencies
#define MATHLIB_STANDALONE

#include <JRmath.h>

#define JWML_WIENER_ERR 1e-10
#define JWML_NAN (0.0 / 0.0)
#define JWML_NEGINF ((-1.0) / 0.0)

double maxOfNum(double, double);
double dwiener(double, double, double, double, double, int);
double pwiener(double, double, double, double, double, int, int);
double qwiener(double, double, double, double, double, int, int);
double rwiener(double, double, double, double, RNG*);

