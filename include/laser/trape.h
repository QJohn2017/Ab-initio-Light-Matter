/* Cathal O Broin - cathal.obroin4 at mail.dcu.ie - 2015
   This work is not developed in affiliation with any organisation.
  
   This file is part of AILM.

   AILM is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   AILM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with AILM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CATHAL_LASER_TRAPE_GUARD
#define CATHAL_LASER_TRAPE_GUARD
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
#include "laser/laser.h"
#include "numeric/type.h"

namespace cathal
{
namespace laser
{
class trape : public pulse<real, real>
{
    protected :
    real Omega, E0, W0, CEP, Main, Ramp, TShift;
    carrier Shape;

    private :

    public :
    trape(unsigned int Train, real Tau, real Shift, real E0, real W0, real CEP, real ramp, real main, carrier shape = cos) : pulse(Train, Tau, Shift), E0(E0), W0(W0), CEP(CEP), Shape(shape)
    {
        Ramp = ramp * 2.0 * pi<real>() / W0;
        Main = main * pi<real>() / W0;
        TShift = Ramp+Main;
        Duration = 2.0*(Ramp + Main) + Train*Tau + Shift;
    }

};

class ltrape : public trape
{
    public :
    using trape::trape;
    real PulseDef(real t)
    {
        t -= TShift;
        return E0 * (fabs(t) <= Main ? 1.0 : (fabs(t) <= Ramp+Main ? 1.0-(fabs(t)-Main)/Ramp : 0.0)) 
                * Shape(W0*t - CEP*pi<real>());
    }
};
class ctrape : public trape
{
    public :
    using trape::trape;
    real PulseDef(real t)
    {
        t -= TShift;
        return E0 * (fabs(t) <= Main ? 1.0 : (fabs(t) <= Ramp+Main ? 0.5*(1.0-cos(M_PI*(fabs(t)-Main-Ramp)/Ramp)) : 0.0)) 
                * Shape(W0*t - CEP*pi<real>());
    }
};
}
}
#endif

