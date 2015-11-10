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
#ifndef CATHAL_LASER_GUARD
#define CATHAL_LASER_GUARD

//Move to math library
#include <boost/math/constants/constants.hpp>
#include "numeric/integrate.h"
#include "numeric/type.h"

using boost::math::constants::pi;

namespace cathal
{
namespace laser
{

typedef real (* carrier)(real);
class pulse
{
    protected :
    unsigned int Train;
    real Tau, Shift, Duration;

    private :
    virtual pulsetype PulseDef(real t)
    {
        return real(0.0);
    }

    public :
    pulse(unsigned int train, real tau, real shift) : Train(train), Tau(tau), Shift(shift) { }
    pulsetype E(real t)
    {
        pulsetype EVal = real(0.0);
        for (unsigned int i = 0; i < Train; i++)
            EVal += PulseDef(t - i*Tau - Shift);
        return EVal;
    }

    friend class field;
};

class field
{
    std::vector<pulse *> Pulses;
    vecpot Vec;
    real Time;
    quadrature::gauss * Gauss;

    public :
    real End;
    pulsetype E(real t);
    vecpot A(real t0, real t1);
    field(real n);
    void AddPulse(pulse * p);
};
pulsetype EVal(real t, field * Pulse)
{
    return Pulse->E(t);
}

pulsetype field::E(real t)
{
    real Val = real(0.0);            
    //for (auto p = Pulses.begin(); p != Pulses.end(); p++)
    for (auto &p : Pulses)
        Val += p->E(t);
}

vecpot field::A(real t0, real t1)
{
    Vec += Gauss->Quad(t0, t1, EVal, this);
}

field::field(real n)
{
    End = real(0.0);
    Vec = vecpot(0.0);
    Gauss = new quadrature::gauss(n);
}

void field::AddPulse(pulse * p)
{
    Pulses.push_back(p);
    if (End < p->Duration)
        End = p->Duration;
}
}
}
#endif

