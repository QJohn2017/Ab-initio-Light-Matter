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
#ifndef CATHAL_LASER_GAUSS_GUARD
#define CATHAL_LASER_GAUSS_GUARD
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
#include "laser/laser.h"
#include "numeric/type.h"

namespace cathal
{
namespace laser
{
class gauss : public pulse<real, real>
{
    private :
    real A, B, C, D, Z0, Zd;
    real E0, W0, CEP;
    real GShift;
    real Duration;
    carrier Shape;
    public :

    gauss(unsigned int Train, real Tau, real Shift, real E0, real W0, real CEP, real FWHM, real d, real Length, carrier Shape = cos) : pulse(Train, Tau, Shift), E0(E0), W0(W0), CEP(CEP), Shape(Shape)
    {
        Z0 = FWHM / (2.0 * sqrt(2.0*log(2.0)));
        Zd = pow(Z0, 4.0) + pow(d, 2.0);
        A = E0 * Z0/pow(Zd, 1.0/4.0);
        B = -pow(Z0, 2.0) / (2.0 * Zd);

        C = -d/(2.0*Zd);
        D = 0.5 * atan(d / pow(Z0, 2.0));

        real TauChirp = 2.0 * sqrt(2.0 * log(2.0) * Zd / (Z0*Z0));
        GShift = Length * TauChirp / 2.0;
        Duration = Length * TauChirp + Shift + Tau*Train;
    }
    real End(void)
    {
        return Duration;
    }
    real PulseDef(real t)
    {
        t -= GShift;
        return A  * exp(B * t*t) * Shape(W0*t - CEP*pi<real>() - (C*t*t + D));
    }
};
}
}
#endif

