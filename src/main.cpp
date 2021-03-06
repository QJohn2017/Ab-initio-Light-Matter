/**
 *
 *  NOTE, this main is not intended to be "clean". It is a test to compile while
 *  developing the libraries. It implements a Hydrogen eigenvalue (and eigenfunction) calculator.
 *
 */

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
#include <iostream>
#include <array>
#include <algorithm>
#include <vector>
#include <utility>
#include <libconfig.h++>
#include <iomanip> //for setprecision

#include "la/array.h"
#include "la/vec.h"
#include "util/io.h"
#include "numeric/integrate.h"
#include "numeric/type.h"
#include "la/krylov.h"
#include "numeric/splines.h"
#include "deprecated.h"
using namespace cathal;
void Resize(std::vector<std::vector<real> > & SplineOverlap, size_t NumKnots, size_t k)
{
    SplineOverlap.resize(NumKnots);
    for (auto & it : SplineOverlap)
        it.resize(2*k-1);
}

namespace cathal
{
    template <class T>
    void ioln(T Out)
    {
        std::cout << Out << std::endl;
    }
}

template <class T>
void MatTest(la::band<T> H, la::band<T> SMat)
{
    ioln("H Eigen");
    SymBandEigen(H);

    la::slice<real> C(H.Row());
    for (auto & t : C)
        t = 5.0;

    Print("CVEC", C);

    ioln("D = H * C;\n");
    la::slice<real> D = H * C;

    Print("DVEC", D);

    ioln("SMat Eigen");
    SymBandEigen(SMat);
}

namespace cathal
{
namespace quant
{
namespace abinitio
{
template <class T>
la::band<T> HOverlapMatrix(quadrature::gauss<T, T> & Gauss, int l, size_t k, std::vector<T> & Knots, size_t IgnoreStart = 1, size_t IgnoreEnd = 1)
{
    size_t No = Knots.size() - k;
    la::band<T> DivX2(No, k), DivX(No, k), H(No, k), Prod(No, k);

    spline::SymmOverlap(Gauss, Knots, k, DivX2, [](real x) {return (x ? 1.0 / (x*x) : 0.0);}, spline::BSpline<real>, spline::BSpline<real>);
    spline::SymmOverlap(Gauss, Knots, k, DivX, [](real x) {return (x ? 1.0 / x : 0.0);}, spline::BSpline<real>, spline::BSpline<real>);
    spline::SymmOverlap(Gauss, Knots, k, H, [&Knots, k](int i, int j, real x) { return spline::DBSpline(k, i, x, Knots) * spline::DBSpline(k, j, x, Knots);});

    for (size_t i = 0; i < Prod.NumElem(); i++)
        Prod(i) = 0.5 * H(i) + 0.5*l*(l+1)*DivX2(i) - DivX(i);

    if (IgnoreStart || IgnoreEnd)
    { 
        la::band<real> Final = Shrink(Prod, IgnoreStart, IgnoreEnd);
        return Final;
    }
    else return Prod;
}
template <class T>
la::band<T> OverlapMatrix(quadrature::gauss<T, T> & Gauss, size_t k, std::vector<T> & Knots, size_t IgnoreStart = 1, size_t IgnoreEnd = 1)
{
    la::band<real> S(Knots.size() - k, k);
    spline::SymmOverlap(Gauss, Knots, k, S, [](real x){return 1.0;}, spline::BSpline<real>, spline::BSpline<real>);
    if (IgnoreStart || IgnoreEnd)
    {
        la::band<real> Final = Shrink(S, IgnoreStart, IgnoreEnd);
        return Final;
    }
    else return S;
}
// std::pair<vec<real>, la:sqrarray<real> >
void Hydrogen(std::vector<real> Knots, int k, int l, size_t NumKrylov, size_t IgnoreStart = 1, size_t IgnoreEnd = 1)
{
    size_t GaussOrder = k;
    quadrature::gauss<real, real> Gauss(GaussOrder);

    la::band<real> H = HOverlapMatrix(Gauss, l, k, Knots, IgnoreStart, IgnoreEnd);
    std::cout << "Hydrogen built\n";
    std::cout << "H.dim = ("<< H.Row() << ", " << H.Column() << ");\n";

    la::band<real> S = OverlapMatrix(Gauss, k, Knots, IgnoreStart, IgnoreEnd);

    la::sqrarray<real> sqrH(1);
    sqrH.AddBlock(0, 0, &H);

    std::cout << std::setprecision(15);
    std::vector<real> Eigen = GenSymBandEigenvalues(H, S);
    ioln("Eigenvalues");

    if (Eigen.size() > 10)
        Eigen.resize(10);
    io::Print(Eigen);

    la::arnoldi<real> Kry(sqrH, NumKrylov);
    ioln("Kry");
    std::vector<real> KryValues = Kry.Eigenvalues();
    io::Print(KryValues);
}
}
}
}
int main(int argc, char ** argv)
{
    std::cout << "SIZES " << sizeof(double) << " " << sizeof(long double) << std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////                    Initialisation                                       ///////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////REPLACE THESE/////////////////
    std::string CFile(argc > 1 ? argv[1] : "settings.cfg");
    size_t NumKrylov = argc > 2 ? atoi(argv[2]) : 5;
///////////////////////////////////////////////////////////////

    std::vector<real> Knots;
    libconfig::Config Conf;
    int k = 3;

    Conf.readFile(CFile.c_str());                           //ParseException, FileIOException
    libconfig::Setting & KnotSet = Conf.lookup("Knots");    //SettingNotFoundException
    libconfig::Setting & kSet = Conf.lookup("k");           //SettingNotFoundException
    std::cout << "Order.\n";
    k = kSet;
    std::cout << "Order = " << k << std::endl;

    Knots.resize(KnotSet.getLength() + 2*k-2);
    for (int i = 0; i < KnotSet.getLength(); i++)
    {
        double Val = KnotSet[i];
        Knots[k-1+i] = real(Val);
    }
    for (int i = 0; i < k-1; i++)
    {
        Knots[i] = Knots[k-1];
        Knots[k + KnotSet.getLength()-1 + i] = Knots[k-1 + KnotSet.getLength()-1];
    }
    io::Print(Knots);
    ioln("\nSettings loaded");

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////                                Main                                     ///////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
    int l = 0;
    quant::abinitio::Hydrogen(Knots, k, l, NumKrylov);

    return 0;
}

