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
#pragma once
#ifndef CATHAL_ERROR_GUARD
#define CATHAL_ERROR_GUARD

namespace cathal
{
typedef enum
{
    UNKNOWN,
    SIZE_MISMATCH,
    DIM_MISMATCH,
    BLOCK_MISMATCH,
    OUT_OF_BOUNDS,
    SELF_ASSIGNMENT
} ErrorCode;


}

#endif

