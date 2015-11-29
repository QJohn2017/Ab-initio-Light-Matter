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
#ifndef CATHAL_KRYLOV_GUARD
#define CATHAL_KRYLOV_GUARD
#include <vector>
#include "la/array.h"
#include "la/vec.h"
#include "deprecated.h"
#include "util/io.h"
namespace cathal
{
namespace la
{
template<class T>
class arnoldi
{
    private :
    bool Defer;
    T Residual = 0.0;
    fullblock<T> HSmall;
    std::vector<vec<T> > QDag;

    public :
    arnoldi(sqrarray<T> & H, size_t N, bool Def = false) : Defer(Def), HSmall(N, N)
    {
        QDag.resize(N+1);

        if (N > H.RowElem())
        {
            DP();
            std::cout << N << " " << H.RowElem() << std::endl;
            throw(UNKNOWN); //Not enough eigenvalues
        }

        for (size_t i = 0; i < N+1; i++)
            QDag[i].Resize(H.Index());

        QDag[0].Set(1.0);
        if (!Defer)
            Residual = Krylov(H);
    }

    std::vector<T> Eigenvalues()
    {
        //This is a temporary place holder for a better method
        la::band<T> TriHSmall = Convert(HSmall, 3);
        return SymBandEigenvalues(TriHSmall);
    }

    T Krylov(sqrarray<T> & H)
    {
        int N = HSmall.Row();
        Normalise(QDag[0]); //Run a normalisation routine

        T Val = 0.0;
        for (int i = 0; i < N; i++)
        {
            vec<T> & P = QDag[i];
            vec<T> & C = QDag[i+1];
            C = slice<T>(H * P); //TODO: Make this cast not required.
            for (int j = 0; j < i+1; j++)
            {
                vec<T> & jth = QDag[j];
                HSmall(j, i) = Dot(jth, C);
                C -= jth * HSmall(j, i);
            }
            Val = Normalise(C); //Run a normalisation routine
            if (i < N-1)
                HSmall(i+1, i) = Val; //Run a normalisation routine
        }
        return Val;
    }
};
}
}
#endif