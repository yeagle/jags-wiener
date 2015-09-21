/*
 *  Copyright (C) 2012 Dominik Wabersich <dominik.wabersich@gmail.com>
 *  and Joachim Vandekerckhove <joachim@uci.edu>
 *
 *  When using this module, please cite as: 
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
#include <module/Module.h>
#include <distributions/DWiener.h>
#include <distributions/DWieners.h>
#include <functions/DFunction.h>
#include <functions/DLogFunction.h>

using std::vector;

namespace jags {
namespace wiener {

class WIENERModule : public Module {
  public:
    WIENERModule();
    ~WIENERModule();
};

WIENERModule::WIENERModule() : Module("wiener")
{
  DWiener *wienerdist;
  wienerdist = new DWiener();
  //load distributions
  insert(wienerdist);
  insert(new DWieners(wienerdist));
  //load functions
  insert(new DFunction(wienerdist));
  insert(new DLogFunction(wienerdist));
}

WIENERModule::~WIENERModule() 
{
  vector<Function*> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
  vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  }
}

} // namespace wiener
} // namespace jags

jags::wiener::WIENERModule _wiener_module;
