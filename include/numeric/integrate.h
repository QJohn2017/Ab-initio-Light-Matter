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
#ifndef CATHAL_QUADRATURE_GUARD
#define CATHAL_QUADRATURE_GUARD

#include <gsl/gsl_integration.h>
#include "numeric/type.h"

namespace cathal
{
namespace quadrature
{
    template <class T=real, class P=real>
    class gauss
    {
        private :
        size_t n;
        //Type is baked in. No choice since the GSL library doesn't provide choice.
        std::vector<double> Xi, Wi;

        void GSLTable(size_t N)
        {
            n = N;
            gsl_integration_glfixed_table * GLTable = gsl_integration_glfixed_table_alloc(n);
            Xi.resize(n);
            Wi.resize(n);

            for (size_t i = 0; i < GLTable->n; i++)
                gsl_integration_glfixed_point(-1.0, 1.0, i, &Xi[i], &Wi[i], GLTable);
            gsl_integration_glfixed_table_free(GLTable);
        }
        
        public :

        gauss(std::vector<std::pair<P, P>> Val) : n(Val.size()), Xi(Val.size()), Wi(Val.size())
        {
            for (size_t i = 0; i < Val.size(); i++)
            {
                Xi[i] = Val[i].first;
                Wi[i] = Val[i].second;
            }
        }

        gauss(size_t N)
        {
            GSLTable(N);
        }

        T Quad(P a, P b, std::function<T(P)> Cast)
        {
            T Val = T(0);
            for (size_t i = 0; i < n; i++)
                Val += Wi[i] * Cast((b-a)/T(2) * Xi[i] + (b + a)/T(2));
            Val *= (b-a)/T(2);
            return Val;
        }
        
        //Recursive calculation.
        T Recursive(T Tol, T Compare, P a, P b, std::function<T(P)> Cast)
        {
            T Val1 = Quad(a, a + (b-a)/P(2), Cast);
            T Val2 = Quad(a + (b-a)/P(2), b, Cast);
            
            if (std::abs(Val1+Val2 - Compare) > Tol)
            {
//                 std::cout << "Compare " << a << " " << b << std::endl;
                T Final = Recursive(Tol, Val1, a, a + (b-a)/P(2), Cast)
                        + Recursive(Tol, Val2, a +  (b-a)/P(2), b, Cast);
                return Final;
            }
            else
                return Val1+Val2;
        }
        T Adaptive(T Tol, P a, P b, std::function<T(P)> Cast)
        {
            T Test = Quad(a, b, Cast);
            return Recursive(Tol, Test, a, b, Cast);
        }
    };
}
}

#endif

