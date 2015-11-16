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
#include <functional>
//Move to math library
// #include <boost/math/constants/constants.hpp>
#include "numeric/integrate.h"
#include "numeric/type.h"

// using boost::math::constants::pi;

namespace cathal
{
namespace laser
{

typedef std::function<real(real)> carrier;

template <class T, class P>
class pulse
{
    protected :
    unsigned int Train;
    T Tau, Shift, Duration;

    private :
    virtual P PulseDef(T t) = 0;
    public :
    T End()
    {
        return Duration;
    }

    pulse(unsigned int train, T tau, T shift) : Train(train), Tau(tau), Shift(shift) { }
    P E(T t)
    {
        P EVal = T(0);
        for (unsigned int i = 0; i < Train; i++)
            EVal += PulseDef(t - i*Tau - Shift);
        return EVal;
    }
};

template <class T, class P>
class field
{
    std::vector<pulse<T, P> *> Pulses;
    P Vec = P(0);
    T Time;
    quadrature::gauss Gauss; //Gaussian quadrature

    public :
    T End = T(0);

    field(T n) : Gauss(n)
    {
    }

    pulsetype E(T t)
    {
        T Val = T(0);
        for (auto p : Pulses)
            Val += p->E(t);
        return Val;
    }

    P A(T t0, T t1)
    {
        Vec += Gauss.Quad<P, T>(t0, t1, [=] (T t) -> P {  return E(t); });
        return Vec;
    }

    void AddPulse(pulse<T, P> * p)
    {
        Pulses.push_back(p);
        if (End < p->End())
            End = p->End();
    }
};
}
typedef laser::field<real, real> rfield;
}
#endif

