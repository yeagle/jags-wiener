/*  This class contains functions for the wiener (diffusion) model */

/*
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
 */
#ifndef DWIENER_H_
#define DWIENER_H_

#include <config.h>
#include <distribution/ScalarDist.h>
#include <JRmath.h>

#define JWML_WIENER_ERR 1e-10
#define JWML_NAN (0.0 / 0.0)
#define JWML_NEGINF ((-1.0) / 0.0)

namespace wiener {

/**
 * @short Diffusion Model distribution
 * <pre>
 * x ~ dwiener(v, a, w)
 * f(x | v, a, w) = ... if x slow
 *                = ... if x fast
 * </pre>
 * Parameters:
 * v - drift rate
 * a - boundary separation
 * w - relative starting point
 */
class DWiener : public ScalarDist 
{
  public:
    DWiener();

    double logDensity(double x, PDFType type,
		      std::vector<double const *> const &parameters,
		      double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
	   	double const *lower, double const *upper,
	   	RNG *rng) const;
    double calPlower(double, std::vector<double const *> const &) const;
    double calPupper(double, std::vector<double const *> const &) const;
    double typicalValue(std::vector<double const *> const &parameters,
	   	double const *lower, double const *upper) const;
    /**
     * Checks that:
     * a > 0
     * w in intervall [0,1] 
     * terr > 0
     */
    bool checkParameterValue(std::vector<double const *> const &parameters) const;

    double d(double x, PDFType type,
      std::vector<double const *> const &parameters,
      bool give_log) const;
    double p(double q, std::vector<double const *> const &parameters, bool lower,
	    bool give_log) const;
    double q(double p, std::vector<double const *> const &parameters, bool lower,
	    bool log_p) const;
    double r(std::vector<double const *> const &parameters, RNG *rng) const;
};

}

#endif /* DWIENER_H_ */
