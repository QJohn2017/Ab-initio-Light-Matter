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
#include <vector>
#include "numeric/type.h"
namespace cathal
{
namespace numeric
{
real BSpline(size_t k, size_t i, real x, std::vector<real> & Knots)
{
    if (k == 1)
        return real(x >= Knots[i] && x < Knots[i+1]);
    if (i+k >= Knots.size())
        return 0.0;

    real A = Knots[i+k-1] - Knots[i];
    real B = Knots[i+k] - Knots[i+1];

    real BS1 = BSpline(k-1, i, x, Knots);
    real BS2 = BSpline(k-1, i+1, x, Knots);

    return  (A ? (x - Knots[i])/A * BS1 : 0.0) + (B ? (Knots[i+k] - x)/B * BS2 : 0.0);
}

real DBSpline(size_t k, size_t i, real x, std::vector<real> & Knots)
{
    if (i+k >= Knots.size())
        return 0.0;

    real A = Knots[i+k-1] - Knots[i];
    real B = Knots[i+k] - Knots[i+1];

    real BS1 = BSpline(k-1, i, x, Knots);
    real BS2 = BSpline(k-1, i+1, x, Knots);

    return (real)(k-1)*((A ? BS1/A : 0.0) - (B ? BS2/B : 0.0));
}
}
}

