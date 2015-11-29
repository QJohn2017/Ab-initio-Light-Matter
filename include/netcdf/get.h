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
/*
* SPECIAL NOTE:
* Requires netcdf-cxx package (archlinux) or equivalent. Uses the (Appel) NetCDF API.
*/
#ifndef NETCDF_GET_GUARD
#define NETCDF_GET_GUARD
#include <ncFile.h>
#include <ncVar.h>
#include <ncDim.h>
//#include <algorithm>
#include <vector>
#include "la/array.h"

namespace cathal
{
namespace nc
{

template <class T>
std::vector<T> GetVector(const std::string & FileName, const std::string & Name)
{
    netCDF::NcFile File(FileName, netCDF::NcFile::read);

    netCDF::NcVar Var = File.getVar(Name);
    std::vector<netCDF::NcDim> Dim = Var.getDims();

    assert(Dim.size() == 1);
    std::vector<T> Data(Dim[0].getSize());
//     Data.resize(Dim[0].getSize());
    Var.getVar(&Data[0]);
    return Data;
}

template <class T>
la::fullblock<T> GetBlock(const std::string & FileName, const std::string & Name)
{
    netCDF::NcFile File(FileName, netCDF::NcFile::read);

    netCDF::NcVar Var = File.getVar(Name);
    std::vector<netCDF::NcDim> Dim = Var.getDims();

    assert(Dim.size() == 2);
    la::fullblock<T> Data(Dim[0].getSize(), Dim[1].getSize());
//     Data.resize(Dim[0].getSize(), Dim[1].getSize());

    Var.getVar(&Data(0));
    return Data;
}
}
}
#endif

