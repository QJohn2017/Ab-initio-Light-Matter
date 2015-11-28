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
#ifndef CATHAL_SLICE_GUARD
#define CATHAL_SLICE_GUARD
#include <vector>
#include <functional>
#include <algorithm>
#include "util/io.h"
namespace cathal
{
namespace la
{
/*
    Slice Requirements:
    A slice denotes a sequence of contiguous memory over which an iterator can iterate from start to finish. Let us use S to represent a particular slice object.
    Converting a la::vec to a la::slice should be well defined.

    Constructors:
    copy constructor is default copy constructor
    slicify la::vec (possible?)
    slicify std::vector
    slicify size (Allocate memory)
    make slice from two iterators

*   Member functions:
    Size(): Number of elements
    SetPair(): To initialise after the fact, accept la::vec, std::vector and iterators

    Also the following operator overloads should be defined:
*   operator* but not *=
    la::block multiplied by S. S multiplied by la::block. (Smoothly treat it as both row/column vector where appropriate)
    S<T> * T
    The la:slice la::vec dot product (and reverse) should resolve to the above dot product.
    The la:slice vector dot product (and reverse) should resolve to the above dot product.
*   operator+ and +=
    addition of another slice
*   operator- and -=
    subtraction of a slice.
*   operator=
    slice, vector, vec should all be supported
*   operator[size_t i]
    return ith element
*   operator/ should not be defined.
*/
/*
 * IDEA: Define slice which has a notion of row vector, column vector
 *
 */
template<class T>
class slice
{
    using iter = typename std::vector<T>::iterator;
    std::vector<T> Mem;
    std::pair<iter, iter> Pair;
    public :
    slice(size_t Num) : Mem(Num)
    {
        SetPair(Mem);
    }
    slice(std::vector<T> & In) : Pair(In.begin(), In.end())
    {
    }
    slice(iter Start, iter End) : Pair(Start, End)
    {
    }
    slice()
    {
    }
 //Initialise after the fact
    void SetPair(iter first, iter second)
    {
        Pair = std::make_pair(first, second);
    }
    void SetPair(std::vector<T> & vIn)
    {
        SetPair(vIn.begin(), vIn.end());
    }

    T & operator[](size_t i)
    {
        if (Pair.first + i < Pair.second)
            return *(Pair.first + i);
        else
            throw(OUT_OF_BOUNDS);
    }

    size_t Size(void)
    {
        return Pair.second - Pair.first;
    }

    iter begin()
    {
        return Pair.first;
    }
    iter end()
    {
        return Pair.second;
    }

    //Operations with two operands
    // A = B Op C
    slice<T> BinOp(slice<T> B, std::function<T(T, T)> Op)
    {
        if (Size() != B.Size())
            throw(SIZE_MISMATCH);

        slice<T> A(Size());
        for (size_t i = 0; i < Size(); i++)
            A[i] = Op(this->operator[](i), B[i]);
        return A;
    }

    //Operations with one operands
    // A = Op C
    slice<T> Unary(std::function<T(T)> Op)
    {
        slice<T> A(Size());
        for (size_t i = 0; i < Size(); i++)
            A[i] = Op(this->operator[](i));
        return A;
    }

//  Binary operator overloads
    slice<T> operator+(slice<T> In)
    {
        return BinOp(In, std::plus<T>());
    }
    slice<T> operator-(slice<T> In)
    {
        return BinOp(In, std::minus<T>());
    }
    slice<T> operator*(T Scal)
    {
        return Unary([Scal](T x) { return Scal * x; });
    }
    slice<T> operator+(T Scal)
    {
        return Unary([Scal](T x) { return Scal * x; });
    }
    slice<T> operator-(T Scal)
    {
        return Unary([Scal](T x) { return Scal - x; });
    }

//  Assignment overload
    slice<T> & operator=(slice<T> In)
    {
        if (In.Size() != Size())
            throw(SIZE_MISMATCH);

        std::copy(In.begin(), In.end(), begin());
        return *this;
    }
    slice<T> & operator=(T Scal)
    {
        std::fill(begin(), end(), Scal);
        return *this;
    }
//  Increment operator overloads
    void operator+=(slice<T> In)
    {
        *this = BinOp(In, std::plus<T>());
    }
    void operator-=(slice<T> In)
    {
        *this = BinOp(In, std::minus<T>());
    }
};
template <class T>
slice<T> operator*(slice<T> A, T Scal)
{
    //TODO: Use std::inner_product ?
    return Scal*slice<T>(A);
}
template <class T>
T Dot(slice<T> A, slice<T> B)
{
    if (A.Size() != B.Size())
        throw(SIZE_MISMATCH);
    T Sum = 0.0;
    for (size_t i = 0; i < A.Size(); i++)
        Sum += std::conj(A[i]) * B[i];
    return Sum;
}
template <class T>
T Normalise(la::slice<T> C)
{
    using iter = typename std::vector<T>::iterator;
    T Sum = 0.0;
    for (iter it = C.begin(); it != C.end(); it++)
        Sum += std::norm(*it);
    Sum = std::sqrt(Sum);
    for (iter it = C.begin(); it != C.end(); it++)
        *it /= Sum;
    return Sum;
}
#include <iostream>
#include <string>
template <class T>
void Print(std::string Name, slice<T> Slice)
{
    std::cout << "std::vector<" << typeid(T).name() << " > " << Name << " = {\n";
        for (size_t j = 0; j < Slice.Size(); j++)
            std::cout << Slice[j] << " (" << j << ")\n";
    std::cout << "};\n";
}
}
}
#endif