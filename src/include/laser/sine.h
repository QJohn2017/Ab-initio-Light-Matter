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

#include "laser/laser.h"
#include "numeric/type.h"

namespace cathal
{
namespace laser
{
class sine : public pulse
{
    protected :
    real Omega, E0, W0, CEP;
    unsigned int Cycles;
    carrier Shape;

    private :

    public :
    sine(unsigned int train, real tau, real shift, real cep, real w0, real e0, unsigned int cycles, carrier shape) : pulse(train, tau, shift),
                CEP(cep), W0(w0), E0(e0), Cycles(cycles), Shape(shape)
    {
        Omega = W0 / (real(2.0) * real(Cycles));
        Duration = pi<real>() / Omega + Train*Tau + Shift;
    }

    pulsetype PulseDef(real t)
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

