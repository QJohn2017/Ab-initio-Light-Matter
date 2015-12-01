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
#ifndef CATHAL_BLOCK_ARRAY_GUARD
#define CATHAL_BLOCK_ARRAY_GUARD
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "util/error.h"
#include "numeric/type.h"
#include "la/slice.h"
#include "la/vec.h"
#include "util/io.h"
namespace cathal
{
namespace la
{
/*
    Block Requirements: (Abstract class)
    A la::block has two Dimensions, Row and Column.

    resize(i, j) should be defined

    operator overload:
    How a la::block multiplies by a la::slice should be defined as a pure virtual function
*/
template <class T, class U=T>
class block
{
    protected:
    size_t N = 0, M = 0;
    public :
    block(size_t N, size_t M) : N(N), M(M)
    {
    }
    size_t Row()
    {
        return N;
    }
    size_t Column()
    {
        return M;
    }
    
    //Note, child classes should throw errors on incorrect access, where incorrect corresponds to "the return value doesn't exist, or would be zero at all times".
    virtual T & operator()(int i) = 0;
    virtual T & operator()(int i, int j) = 0;
    virtual size_t NumElem() = 0;
    virtual slice<U> operator*(slice<U> A) = 0;
};
template <class T, class U=T>
class fullblock : public block<T, U>
{
    std::vector<T> Data;
    public :
    fullblock(size_t N, size_t M) : block<T, U>(N, M), Data(N*M)
    {
    }
    fullblock(std::initializer_list<T> Init, size_t N, size_t M)
    {
        fullblock(N, M);
        if (N*M < Init.size())
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        std::copy(Init.begin(), Init.end(), Data.begin());
    }
    fullblock(std::vector<T> Init, size_t N, size_t M) : block<T, U>(N, M), Data(N*M)
    {
        fullblock(N, M);
        if (N*M < Init.size())
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        std::copy(Init.begin(), Init.end(), Data.begin());
    }
    void Resize(size_t i, size_t j)
    {
        this->N = i;
        this->M = j;
        Data.resize(this->N*this->M);
    }
    slice<U> operator*(slice<U> B)
    {
        slice<U> A(this->N);
        if (B.Size() != this->M)
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        for (size_t i = 0; i < this->N; i++)
            for (size_t j = 0; j < this->M; j++)
                    A[i] += this->operator()(i, j) * B[j];
        return A;
    }
    //For working out calculation complexity
    size_t NumElem()
    {
        return this->N*this->M;
    }
    T & operator()(int i) // const
    {
        return Data.at(i);
    }
    T & operator()(int i, int j) // const
    {
        return Data.at(i*this->M + j);
    }
};

template <class T, class U=T>
class band : public block<T, U>
{
    int k;
    std::vector<T> Data;
    public :
    //For working out calculation complexity
    size_t NumElem()
    {
        return this->N*(2*k-1);
    }
    band(size_t N, int k) : block<T, U>(N, N), k(k), Data(N*(2*k-1))
    {
    }
    band(std::initializer_list<T> Init, size_t N, int k) : band(N, k)
    {
        if (NumElem() != Init.size())
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        std::copy(Init.begin(), Init.end(), Data.begin());
    }
    /*
    band(std::vector<T> Init, size_t N, size_t M) : block<T, U>(N, M), Data(N*M)
    {
        band(k, N);
        if (N*M != Init.size())
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        std::copy(Init.begin(), Init.end(), Data.begin());
    }*/
    int Order()
    {
        return k;
    }
    void Resize(size_t i, int Newk)
    {
        this->N = this->M = i;
        k = Newk;
        Data.resize(NumElem());
    }
    slice<U> operator*(slice<U> B)
    {
        size_t Sz = this->Row(); //Square matrix, Rows() == Columns()
        slice<U> A(Sz);
        if (B.Size() != Sz)
            throw(SIZE_MISMATCH);
        for (int i = 0; i < int(Sz); i++)
            for (int j = std::max(0, i-k+1); j < std::min(int(Sz), i+k); j++)
                    A[i] += this->operator()(i, j) * B[j];
        return A;
    }

    T & operator()(int i) // const
    {
        return Data.at(i);
    }
    T & operator()(int i, int j) // const
    {
//         return Data.at(size_t(i*(2*k-1) + j - i + k - 1));
        return Data.at(size_t((i*2+1)*(k-1) + j));
    }
};

template <class T, class U=T>
class diag : public block<T, U>
{
    int k;
    std::vector<T> Data;
    public :
    //For working out calculation complexity
    size_t NumElem()
    {
        return this->N;
    }
    diag(std::vector<real> & Init, size_t N) : block<T, U>(N, N)
    {
        Data = Init;
    }
    diag(std::initializer_list<T> Init, size_t N) : diag(N)
    {
        if (NumElem() != Init.size())
        {
            DP();
            throw(SIZE_MISMATCH);
        }
        std::copy(Init.begin(), Init.end(), Data.begin());
    }
    void Resize(size_t i)
    {
        this->N = this->M = i;
        Data.resize(NumElem());
    }
    slice<U> operator*(slice<U> B)
    {
        size_t Sz = this->Row(); //Square matrix, Rows() == Columns()
        slice<U> A(Sz);
        if (B.Size() != Sz)
            throw(SIZE_MISMATCH);
        for (int i = 0; i < int(Sz); i++)
            A[i] += this->operator()(i) * B[i];
        return A;
    }

    T & operator()(int i) // const
    {
        return Data.at(i);
    }
    T & operator()(int i, int j) // const
    {
        if (i == j)
            return Data.at(size_t(i));
        else
            throw (OUT_OF_BOUNDS);  //This could be trivially compensated for, but it's really a sign of an algorithm problem.
    }
};

template <class T>
band<T> Convert(fullblock<T> & B, int k)
{
    band<T> A(B.Row(), k);
    for (size_t i = 0; i < A.Row(); i++)
        for (int j = std::max(0, int(i)-k+1); j < std::min(int(A.Row()), int(i)+k); j++)
            A(i, j) = B(i, j);
    return A;
}

template <class T, class U>
struct BlockElement
{
    size_t i;
    size_t j;
    block<T, U> * Arg;
};

//T is the type of the matrix elements, U is the type of Supported vector multiplication
template <class T, class U = T>
class sqrarray //: public block<T, U> //TODO: Add slice multiply support to fufill block class. This will allow blocks of blocks.
{
    size_t N, M;
    bool mirror = false; //An access to j,i will be converted to i,j
    public :
    std::vector<BlockElement<T, U> > Data;
    std::vector<size_t> Indx;
    sqrarray(size_t N) : N(N), M(N), Indx(N)
    {

    }
    void AddBlock(size_t i, size_t j, block<T, U> * Arg)
    {
        Data.push_back({i, j, Arg});
        Indx.at(i) = Arg->Row(); 
        Indx.at(j) = Arg->Column(); 
    }
    size_t Row()
    {
        return N;
    }
    size_t RowElem()
    {
        size_t Sum = 0;
        for (auto i : Indx)
            Sum += i;
        return Sum;
    }
    std::vector<size_t> & Index()
    {
        return Indx;
    }
    size_t Column()
    {
        return Row();
    }
    size_t ColumnElem()
    {
        return RowElem();
    }
    size_t NumElem()
    {
        size_t Sum = 0;
        //Iterate through blocks
        for (auto & k : Data)
            Sum += k.Arg->NumElem();
        return Sum;
    }
    vec<U> operator*(vec<U> & B)
    {
        vec<U> A(B.Index());    //Square matrix, so this is ok
        if (B.Blocks() == Column())
        for (auto k : Data)
        {
            slice<U> sA = A.Block(k.i);
            slice<U> sB = B.Block(k.j);
            sA = *k.Arg * sB;
        }
        else
            throw(BLOCK_MISMATCH);
        return A;
    }
};
template <class T>
la::band<T> Shrink(la::band<T> Full, size_t Start, size_t End)
{
//TODO Check Row()-End-Start > 0?
    int k = Full.Order();
    size_t Sz = Full.Row();
    la::band<T> Shrunk(Sz-Start-End, k);
    for (int i = Start; i < int(Sz-End); i++)
        for (int j = std::max(int(Start), i-k+1); j < std::min(int(Sz-End), i+k); j++)
            Shrunk(i-Start, j-Start) = Full(i, j);
    return Shrunk;
}
}
}
#endif
