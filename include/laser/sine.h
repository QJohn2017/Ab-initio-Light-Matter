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
#ifndef CATHAL_LASER_SINE_GUARD
#define CATHAL_LASER_SINE_GUARD
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
#include "laser/laser.h"
#include "numeric/type.h"

namespace cathal
{
namespace laser
{
class sine : public pulse<real, real>
{
    protected :
    real Omega, E0, W0, CEP;
    unsigned int Cycles;
    std::function<real(real)> Shape;

    private :

    public :
    sine(unsigned int Train, real Tau, real Shift, real E0, real W0, real CEP, unsigned int Cycles, carrier Shape = cos) : pulse<real, real>(Train, Tau, Shift), E0(E0), W0(W0), CEP(CEP), Cycles(Cycles), Shape(Shape)
    {
        Omega = W0 / (real(2.0) * real(Cycles));
        Duration = pi<real>() / Omega + Train*Tau + Shift;
    }

    real PulseDef(real t)
    {
        real SineSqr;
        if (t > pi<real>() / Omega || t < real(0.0))
            return real(0.0);

        SineSqr = sin(Omega*t);
        SineSqr *= SineSqr;
        return E0 * SineSqr * Shape(W0*t- CEP*pi<real>());
    }
};
}
}
#endif

