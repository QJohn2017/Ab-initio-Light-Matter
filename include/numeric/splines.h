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
#ifndef SPLINES_INCLUDE_GUARD
#define SPLINES_INCLUDE_GUARD

#include <vector>
#include <functional>
#include "numeric/type.h"
#include "util/array.h"
namespace cathal
{
namespace numeric
{
// typedef real (* QSOLFunc)(real, void *);
typedef std::function<real(real)> QSOLFunc;

// extern int NumPoints;
// extern std::vector<real> Xi;
// extern std::vector<real> Wi;
// extern void InitQuad(int n);
//extern real * Wi;

extern real BSpline(size_t k, size_t i, real x, std::vector<real> & Knots);
extern real DBSpline(size_t k, size_t i, real x, std::vector<real> & Knots);

typedef std::function<real(size_t, size_t, real, std::vector<real> &)> spl;
typedef std::function<real(int, int, real)> QuadFunc;

void BSplineOverlap(quadrature::gauss & Gauss, std::vector<real> & Knots, int k, util::array<real> & SplineOverlap, QuadFunc Fun)
{
    using namespace std::placeholders;
    int Ns = SplineOverlap.dim(0);
#pragma omp parallel for shared(Fun, Gauss, SplineOverlap, Knots) firstprivate(Ns, k) default(none)
    for (int i = 0; i < Ns; i++)
        for (int j = std::max(0, i-k+1); j < std::min(Ns, i+k); j++)
        {
            real Sum = 0.0;
            int Min = std::max(std::min(i, j) - k+1, 0);
            for (int bps = Min; bps < Min+2*k-1; bps++)
                Sum += Gauss.Quad<real, real>(Knots[bps], Knots[bps+1], std::bind(Fun, i, j, _1));
            SplineOverlap(i, j) = Sum;
        }
}

//Same functionality as above except optimised for symmetric problems (when <i|O|j> == <j|O|i>)
//Tested on k=3, 100,000 by 100,000 for speed, ~30% reduction in run time (openmp disabled)
//Tested on k=9, 1,000 by 1,000 for speed, ~30% reduction in run time (openmp disabled)
void SymmBSplineOverlap(quadrature::gauss & Gauss, std::vector<real> & Knots, int k, util::array<real> & SplineOverlap, QuadFunc Fun)
{
    using namespace std::placeholders;
    int Ns = SplineOverlap.dim(0);
    
#pragma omp parallel for shared(Fun, Gauss, SplineOverlap, Knots) firstprivate(Ns, k) default(none)
    for (int i = 0; i < Ns; i++)
    {
        int Min = std::max(i - k+1, 0);
        
        real Sum = 0.0;
        for (int bps = Min; bps < Min+2*k-1; bps++)
            Sum += Gauss.Quad<real, real>(Knots[bps], Knots[bps+1], std::bind(Fun, i, i, _1));
        SplineOverlap(i, i) = Sum;

        for (int j = i+1; j < std::min(Ns, i+k); j++)
        {
            real Sum = 0.0;
            for (int bps = Min; bps < Min+2*k-1; bps++)
                Sum += Gauss.Quad<real, real>(Knots[bps], Knots[bps+1], std::bind(Fun, i, j, _1));
            SplineOverlap(i, j) = SplineOverlap(j, i) = Sum;
        }
    }
}

void SymmBSplineOverlap(quadrature::gauss & Gauss, std::vector<real> & Knots, int k, util::array<real> & SplineOverlap, QSOLFunc Fun, spl Spl1, spl Spl2)
{
    SymmBSplineOverlap(Gauss, Knots, k, SplineOverlap, [&](int i, int j, real x) -> real { return Fun(x) * Spl1(k, i, x, Knots) * Spl2(k, j, x, Knots);});
}
void BSplineOverlap(quadrature::gauss & Gauss, std::vector<real> & Knots, int k, util::array<real> & SplineOverlap, QSOLFunc Fun, spl Spl1, spl Spl2)
{
    BSplineOverlap(Gauss, Knots, k, SplineOverlap, [&](int i, int j, real x) -> real { return Fun(x) * Spl1(k, i, x, Knots) * Spl2(k, j, x, Knots);});
}

/* BASIS NAIVE IMPLEMENTATION */ /*
template<class T>
real EvalSplineCoef(int k, real x, std::vector<real> & Knots, std::vector<T> & SplineCoef)
{
     real Sum = 0.0;
     for (int i = 0; i < SplineCoef.size(); i++)
        Sum += BSpline(k, i, x, Knots) * SplineCoef[i];
    return Sum;
}

// */

template<class T>
real EvalSplineCoef(int k, real x, std::vector<real> & Knots, std::vector<T> & SplineCoef)
{
     real Sum = 0.0;

    size_t j = k-1;
    while ((Knots[j] > x || (Knots[j+1] < x)) && j+1 < SplineCoef.size()) j++;

    for (int i = j-k+1; i < std::min(j+k-1, SplineCoef.size()); i++)
        Sum += BSpline(k, i, x, Knots) * SplineCoef[i];
    return Sum;
}
}
}
#endif

