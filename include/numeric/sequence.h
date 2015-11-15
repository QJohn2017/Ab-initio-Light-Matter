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
#ifndef CATHAL_SEQ_GUARD
#define CATHAL_SEQ_GUARD
#include "numeric/type.h"
namespace cathal
{
    namespace sequences
    {
        class sequence
        {
            protected :
            real Duration;
            public :
            sequence(real end) : Duration(end) { }

            virtual real Next()
            {
                return real(0.0);
            }
            virtual bool End()
            {
                return true;
            }
        };
        class linear : public sequence
        {
            private :
            unsigned int StepNum;
            real dt;
            public :
            linear (real end, unsigned int NumSteps) : sequence(end)
            {
                StepNum = 0;
                dt = Duration / real(NumSteps+1);
            }

            real Next()
            {
                return dt*StepNum++;
            }
            bool End()
            {
                return (dt*StepNum >= Duration);
            }
        };
    }
}
#endif

