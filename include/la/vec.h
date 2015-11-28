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
#ifndef CATHAL_BLOCK_VEC_GUARD
#define CATHAL_BLOCK_VEC_GUARD
#include <vector>
#include "util/error.h"
#include "la/slice.h"
namespace cathal
{
namespace la
{
/* ***************************************************
 *
 *          Linear Algebra Routines
 *                 -Vector-
 * ***************************************************/
/*
    la::vec Requirements:
    A la::vec denotes a structured vector. Let us use V to represent a particular la::vec object.
    Converting a la::vec to a la::slice should be well defined.
    Converting a block of la::vec to a la::slice should be well defined

    Functions:
    The la:vec la::slice dot product (and reverse) should resolve to a la::slice la::slice dot product.
    V * la::vec should resolve to a la::slice la::slice dot product.

    Also the following operator overloads should be defined:
*   operator*
    la::block multiplied by V. V multiplied by la::block. (Smoothly treat it as both row/column vector where appropriate)
    la:vec std::vector (and reverse) should resolve to a la::slice la::slice dot product.
*   operator+
    addition of another la::vec
*   operator-
    subtraction of a la::vec.
*   operator=
    slice, vector, vec should all be supported
*   operator/ should not be defined.
*/
template <class T>
class vec
{
    using iter = typename std::vector<T>::iterator;

    std::vector<T> Mem;
    slice<T> DataS;
    std::vector<slice<T> > Slices;
    std::vector<size_t> Indx;
    std::vector<size_t> StepIn;
    public :
    void Resize(std::vector<size_t> & In) //Assumes nothing about the current status
    {
        size_t Nb = In.size();  //Number of blocks
        StepIn.resize(Nb + 1);
        Indx = In;

        StepIn[0] = 0;
        for (size_t i = 0; i < Nb; i++)
            StepIn[i+1] = StepIn[i] + Indx[i];

        Mem.resize(StepIn.back());
        DataS.SetPair(Mem);
        Slices.resize(Nb);

        iter it = Mem.begin();
        for (size_t i = 0; i < Nb; i++)
        {
            Slices[i].SetPair(it, it + Indx[i]);
            it += Indx[i];
        }
    }
    vec(){}
    
    vec(size_t Sz)   //A vector consisting of only one logical block
    {
        std::vector<size_t> In = {Sz};
        Resize(In);
    }

    vec(std::vector<size_t> & In)
    {
        Resize(In);
    }

    std::vector<size_t> & Index()
    {
        return Indx;
    }

    slice<T> Block()
    {
        return DataS;
    }
    slice<T> Block(size_t i)
    {
        return Slices[i];
    }
    size_t Blocks()
    {
        return Indx.size();
    }
    size_t Size()
    {
        return StepIn.back();
    }
    slice<T> operator[](size_t i)
    {
        return Block(i);
    }
    T & operator()(size_t i)
    {
        return Mem[i];
    }
//operator overloads
    operator slice<T>()
    {
        return Block();
    }
//  Binary operator overloads
    vec<T> operator+(vec<T> & In)
    {
        vec<T> Out(In.Index());
        Out = Block() + slice<T>(In);
        return Out;
    }
    vec<T> operator-(vec<T> & In)
    {
        vec<T> Out(In.Index());
        Out = Block() - slice<T>(In);
        return Out;
    }
    vec<T> operator*(T Scal)
    {
        vec<T> Out(Index());
        Out = Block() * Scal;
        return Out;
    }
    vec<T> operator+(T Scal)
    {
        vec<T> Out(Index());
        Out = Block() + Scal;
        return Out;
    }
    vec<T> operator-(T Scal)
    {
        vec<T> Out(Index());
        Out = Block() - Scal;
        return Out;
    }
//  Assignment overload
    vec<T> & operator=(vec<T> B)
    {
        if (&B == this)
            throw(SELF_ASSIGNMENT); //I don't see why this would ever be desirable. Change when you see a use case.

        std::swap(Mem, B.Mem);
        std::swap(DataS, B.DataS);
        std::swap(Slices, B.Slices);
        std::swap(Indx, B.Indx);
        std::swap(StepIn, B.StepIn);
        return *this;
    }
    vec<T> & operator=(slice<T> B)
    {
        std::copy(B.begin(), B.end(), Block().begin());
        return *this;
    }
    void Set(T Scal)
    {
        Block() = Scal;
    }

    void operator+=(slice<T> In)
    {
        Block() += In;
    }
    void operator-=(slice<T> In)
    {
        Block() -= slice<T>(In);
    }
};
template <class T>
vec<T> operator*(T Scal, vec<T> & A)
{
    return Scal*A;
}

template <class T>
T Dot(vec<T> & A, vec<T> & B)
{
    return Dot(slice<T>(A), slice<T>(B));
}
template <class T>
T Normalise(la::vec<T> & C)
{
    return Normalise(slice<T>(C));
}
#include <iostream>
#include <string>
template <class T>
void Print(std::string Name, vec<T> Vec)
{
    std::cout << "std::vector<" << typeid(T).name() << " > " << Name << " = {\n";
    for (size_t i = 0; i < Vec.Blocks(); i++)
        for (size_t j = 0; j < Vec.Block(i).Size(); j++)
            std::cout << Vec[i][j] << " (" << i << ", " << j << ")\n";
    std::cout << "};\n";
}
}
}
#endif
