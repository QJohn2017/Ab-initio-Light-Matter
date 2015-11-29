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
#include "la/array.h"
namespace cathal
{
namespace spline
{
// typedef real (* QSOLFunc)(real, void *);
typedef std::function<real(real)> QSOLFunc;

// extern real BSpline(size_t k, size_t i, real x, std::vector<real> & Knots);
// extern real DBSpline(size_t k, size_t i, real x, std::vector<real> & Knots);

typedef std::function<real(size_t, size_t, real, std::vector<real> &)> spl;
typedef std::function<real(int, int, real)> QuadFunc;

template <class T>
T BSpline(size_t k, size_t i, real x, std::vector<T> & Knots)
{
    if (k == 1)
        return real(x >= Knots[i] && x < Knots[i+1]);
    if (i+k >= Knots.size())
        return 0.0;

    T A = Knots[i+k-1] - Knots[i];
    T B = Knots[i+k] - Knots[i+1];

    T BS1 = BSpline(k-1, i, x, Knots);
    T BS2 = BSpline(k-1, i+1, x, Knots);

    return  (A ? (x - Knots[i])/A * BS1 : 0.0) + (B ? (Knots[i+k] - x)/B * BS2 : 0.0);
}

template <class T>
T DBSpline(size_t k, size_t i, real x, std::vector<T> & Knots)
{
    if (i+k >= Knots.size())
        return 0.0;

    T A = Knots[i+k-1] - Knots[i];
    T B = Knots[i+k] - Knots[i+1];

    T BS1 = BSpline(k-1, i, x, Knots);
    T BS2 = BSpline(k-1, i+1, x, Knots);

    return (real)(k-1)*((A ? BS1/A : 0.0) - (B ? BS2/B : 0.0));
}

template <class T, class P>
void Overlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QuadFunc Fun)
{
    using namespace std::placeholders;
    int Ns = SplineOverlap.Row();
#pragma omp parallel for shared(Fun, Gauss, SplineOverlap, Knots) firstprivate(Ns, k) default(none)
    for (int i = 0; i < Ns; i++)
        for (int j = std::max(0, i-k+1); j < std::min(Ns, i+k); j++)
        {
            real Sum = 0.0;
            int Min = std::max(std::min(i, j) - k+1, 0);
            for (int bps = Min; bps < Min+2*k-1; bps++)
                Sum += Gauss.Quad(Knots[bps], Knots[bps+1], std::bind(Fun, i, j, _1));
            SplineOverlap(i, j) = Sum;
        }
}

//Same functionality as above except optimised for symmetric problems (when <i|O|j> == <j|O|i>)
//Tested on k=3, 100,000 by 100,000 for speed, ~30% reduction in run time (openmp disabled)
//Tested on k=9, 1,000 by 1,000 for speed, ~30% reduction in run time (openmp disabled)
template <class T, class P>
void SymmOverlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QuadFunc Fun)
{
    using namespace std::placeholders;
    int Ns = SplineOverlap.Row();

#pragma omp parallel for shared(Fun, Gauss, SplineOverlap, Knots) firstprivate(Ns, k) default(none)
    for (int i = 0; i < Ns; i++)
    {
        int Min = std::max(i - k+1, 0);

        real Sum = 0.0;
        for (int bps = Min; bps < Min+2*k-1; bps++)
            Sum += Gauss.Quad(Knots[bps], Knots[bps+1], std::bind(Fun, i, i, _1));
        SplineOverlap(i, i) = Sum;

        for (int j = i+1; j < std::min(Ns, i+k); j++)
        {
            real Sum = 0.0;
            for (int bps = Min; bps < Min+2*k-1; bps++)
                Sum += Gauss.Quad(Knots[bps], Knots[bps+1], std::bind(Fun, i, j, _1));
            SplineOverlap(i, j) = SplineOverlap(j, i) = Sum;
        }
    }
}
template <class T, class P>
void AdaptiveOverlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QuadFunc Fun)
{
    using namespace std::placeholders;
    int Ns = SplineOverlap.Row();
    #pragma omp parallel for shared(std::cout, Fun, Gauss, SplineOverlap, Knots) firstprivate(Ns, k) default(none)
    for (int i = 0; i < Ns; i++)
    {
        for (int j = std::max(0, i-k+1); j < std::min(Ns, i+k); j++)
        {
            int Start = std::max(std::min(i, j) - k+1, 0);
            int End = Start + 2*k-1;
            real Sum = 0.0;

            SplineOverlap(i, j) = Gauss.Adaptive(1e-16, Knots[Start], Knots[End], std::bind(Fun, i, j, _1));
        }
    }
}
template <class T, class P>
void AdaptiveOverlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QSOLFunc Fun, spl Spl1, spl Spl2)
{
    AdaptiveOverlap(Gauss, Knots, k, SplineOverlap, [&](int i, int j, real x) -> real { return Fun(x) * Spl1(k, i, x, Knots) * Spl2(k, j, x, Knots);});
}
template <class T, class P>
void SymmOverlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QSOLFunc Fun, spl Spl1, spl Spl2)
{
    SymmOverlap(Gauss, Knots, k, SplineOverlap, [&](int i, int j, real x) -> real { return Fun(x) * Spl1(k, i, x, Knots) * Spl2(k, j, x, Knots);});
}
template <class T, class P>
void Overlap(quadrature::gauss<T, P> & Gauss, std::vector<P> & Knots, int k, la::block<T> & SplineOverlap, QSOLFunc Fun, spl Spl1, spl Spl2)
{
    Overlap(Gauss, Knots, k, SplineOverlap, [&](int i, int j, real x) -> real { return Fun(x) * Spl1(k, i, x, Knots) * Spl2(k, j, x, Knots);});
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

