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
/*****************************************************************
 *
 *  Unit tests for core functionality
 *
 *****************************************************************/
#include <iostream>
#include <vector>
#include <array>

//This is to prevent type.h changing the type.
#define CATHAL_TYPE_GUARD
typedef double real;

#include "laser/sine.h"
#include "laser/gauss.h"
#include "laser/trape.h"
#include "la/array.h"
#include "la/vec.h"
#include "la/slice.h"
#include "la/krylov.h"
#include "util/io.h"
#include "numeric/sequence.h"
#include "numeric/splines.h"
#include <gtest/gtest.h>


using namespace cathal;
using std::cout;
using std::endl;

//TODO: Add death tests
void PrintLoop(rfield & Laser)
{
    int NumSteps = 10;
    sequences::sequence * seq = new sequences::linear(Laser.End, NumSteps);
    real t = 0.0, Time = 0.0;
    std::cout << "hexfloat:\n" << std::hexfloat;
    cout << "std::vector<std::array<real, 3>> Test {";

    bool once = true;
    while (!seq->End())
    {
        t = seq->Next();
        real E = Laser.E(t);
        real A = Laser.A(Time, t);

        if (!once)
            cout << ",";
        else
            once = false;

        cout << "{" << t << ", " << E << ", " << A << "}";

        Time = t;
    }
    cout << "};\n";
    cout << std::defaultfloat;
}

std::vector<std::array<real, 3>> Loop(rfield & Laser)
{
    std::vector<std::array<real, 3>> Values;
    int NumSteps = 10;
    sequences::sequence * seq = new sequences::linear(Laser.End, NumSteps);
    real t = 0.0, Time = 0.0;

    int i = 0;
    while (!seq->End())
    {
        t = seq->Next();
        real E = Laser.E(t);
        real A = Laser.A(Time, t);

        Values.push_back({t, E, A});
        Time = t;
        i++;
    }
    return Values;
}
void Compare(std::vector<std::array<real, 3>> x, std::vector<std::array<real, 3>> y)
{
    EXPECT_EQ(x.size(), y.size()) << "Number of time steps do not match\n";
    for (size_t i = 0; i < x.size(); i++)
    {
        ASSERT_DOUBLE_EQ(x[i][0], y[i][0]) << "Time Differs: " << i << std::endl;
        ASSERT_DOUBLE_EQ(x[i][1], y[i][1]) << "E(t) Differs: " << i << std::endl;
        ASSERT_DOUBLE_EQ(x[i][2], y[i][2]) << "A(t) Differs: " << i << std::endl;
    }
}
void Compare(std::vector<real> x, std::vector<real> y)
{
    EXPECT_EQ(x.size(), y.size());
    for (size_t i = 0; i < x.size(); i++)
        ASSERT_DOUBLE_EQ(x[i], y[i])  << "Loop: " << i << std::endl;
}

void Compare(la::band<real> x, std::vector<real> y)
{
    EXPECT_EQ(x.NumElem(), y.size());
    for (size_t i = 0; i < x.NumElem(); i++)
        ASSERT_DOUBLE_EQ(x(i), y[i])  << "Loop: " << i << std::endl;
}

TEST(LaserPulse, HandlesSine)
{
    SCOPED_TRACE("Sine test\n");

    real cep = 0.2;
    real w0 = 1.0;
    real e0 = 0.2;
    unsigned int train = 1;
    unsigned int cycles = 10;
    laser::carrier shape = sin;

    unsigned int GaussN = 9;
    rfield Laser(GaussN);

    Laser.AddPulse(new laser::sine(train, 0.0, 0.0, cep, w0, e0, cycles, shape));
    std::vector<std::array<real, 3>> Values = Loop(Laser);

    std::vector<std::array<real, 3>> Test {{0x0p+0, 0x0p+0, 0x0p+0},{0x1.6d91306c99d5cp+2, -0x1.e4bcc7ea6644bp-7, -0x1.7427a1631f65fp-7},{0x1.6d91306c99d5cp+3, -0x1.d55a8d50366ap-5, 0x1.be959c7c20758p-10},{0x1.122ce45173605p+4, -0x1.4f8a543722de2p-4, 0x1.25c94ee2094bap-4},{0x1.6d91306c99d5cp+4, -0x1.330d9419896d3p-5, 0x1.463668c64131p-3},{0x1.c8f57c87c04b3p+4, 0x1.0db74e5ef42aep-4, 0x1.7c1ed703bec8cp-3},{0x1.122ce45173605p+5, 0x1.3dcb2db878c6bp-3, 0x1.dff5138fbd1ffp-4},{0x1.3fdf0a5f069bp+5, 0x1.51ac1841a1178p-3, 0x1.7bc8f0222ca1p-8},{0x1.6d91306c99d5cp+5, 0x1.9dd1b93832e0bp-4, -0x1.02311a0c8535ap-4},{0x1.9b43567a2d108p+5, 0x1.d6f954b044721p-6, -0x1.c9b50cb2a4a88p-5},{0x1.c8f57c87c04b3p+5, -0x1.db24ca1e8a6d2p-11, -0x1.fc183977816a4p-7}};

    Compare(Values, Test);
}
TEST(LaserPulse, HandlesGauss)
{
    SCOPED_TRACE("Gauss test\n");
    real cep = 0.0;
    real w0 = 1.0;
    real e0 = 0.2;
    real tau = 0.0;
    real shift = 10.0;
    unsigned int train = 1;
    laser::carrier shape = sin;
    unsigned int GaussN = 9;

    rfield Laser(GaussN);
    Laser.AddPulse(new laser::gauss(train, tau, shift, e0, w0, cep, 10.0, 0.0, 5.0));
//      PrintLoop(Laser);
    std::vector<std::array<real, 3>> Values = Loop(Laser);

    std::vector<std::array<real, 3>> Test {{0x0p+0, -0x1.7226feea6dd5ap-52, 0x0p+0},{0x1.5d1745d1745d1p+2, -0x1.001828467497cp-39, 0x1.a2b4dbea6fde9p-41},{0x1.5d1745d1745d1p+3, 0x1.6431bc2691813p-27, 0x1.80274918dd99cp-27},{0x1.05d1745d1745dp+4, 0x1.af47090a1b8e9p-17, 0x1.082bf8402990cp-17},{0x1.5d1745d1745d1p+4, 0x1.5a23bd40c8c44p-10, 0x1.1cb89dfb2632cp-14},{0x1.b45d1745d1745p+4, 0x1.3c5922461a20bp-8, -0x1.f2e54116a051ap-6},{0x1.05d1745d1745dp+5, -0x1.ca5fbec696597p-4, -0x1.3dd31617e0c19p-3},{0x1.31745d1745d17p+5, -0x1.351a784f59dcdp-3, 0x1.9ffe7428e4cdp-6},{0x1.5d1745d1745d1p+5, -0x1.24155602ab919p-6, 0x1.7318ac4e1371ep-6},{0x1.88ba2e8ba2e8bp+5, 0x1.3b7240dbda55ap-15, 0x1.8004ba6ea5cep-11},{0x1.b45d1745d1745p+5, 0x1.02b9ad4e4d7fep-18, 0x1.0e5808300682ep-12},{0x1.dffffffffffffp+5, 0x1.95ff25e10d10fp-28, 0x1.0ee3948570436p-12}};
    Compare(Values, Test);
}

TEST(LaserPulse, HandlesTrape)
{
    SCOPED_TRACE("Trapezoidal test\n");
    real cep = 0.0;
    real w0 = 1.0;
    real e0 = 0.2;
    real tau = 0.0;
    real shift = 10.0;
    unsigned int train = 1;
    laser::carrier shape = sin;
    unsigned int GaussN = 9;

    rfield Laser(GaussN);

    Laser.AddPulse(new laser::ctrape(train, tau, 2.0*shift, e0, w0, cep, 2.0, 10.0, shape));
    std::vector<std::array<real, 3>> Values = Loop(Laser);

    std::vector<std::array<real, 3>> Test {{0x0p+0, 0x0p+0, 0x0p+0},{0x1.3a142d88879c9p+3, 0x0p+0, 0x0p+0},{0x1.3a142d88879c9p+4, 0x0p+0, 0x0p+0},{0x1.d71e444ccb6aep+4, -0x1.c2e79a36b2739p-9, 0x1.598a7c84e6e6ap-3},{0x1.3a142d88879c9p+5, 0x1.46c377ddcdabdp-4, -0x1.92a05d90c1128p-3},{0x1.889938eaa983bp+5, -0x1.25f76ba8cafbfp-3, 0x1.023512a3f1166p-3},{0x1.d71e444ccb6aep+5, 0x1.7c5c8244b1613p-3, -0x1.66003a03aff64p-4},{0x1.12d1a7d776a9p+6, -0x1.99944c7104eb8p-3, -0x1.f2413614b0e18p-7},{0x1.3a142d88879c9p+6, 0x1.793a6855ee701p-3, 0x1.09261c8c052f1p-4},{0x1.6156b33998902p+6, -0x1.202bcabf2c337p-3, -0x1.3e1b0cde51f66p-3},{0x1.889938eaa983bp+6, 0x1.142a225dd0acdp-4, 0x1.3f12f85d7a258p-3},{0x1.afdbbe9bba774p+6, 0x0p+0, 0x1.a8b42f186f8p-13}};
     Compare(Values, Test);
}

//Hamiltonian post-overlap matrix calculation for the test case.
la::band<real> Ham({0, 0, 0x1.83e4e8f93a3e4p+2, -0x1.1104104104102p+0, -0x1.ad8b8362e0d8bp-1,
0, -0x1.1104104104103p+0, 0x1.26b560d826b09p+2, -0x1.c73b23f5651bcp+0, -0x1.ac615d2ace7a9p-1,
-0x1.ad8b8362e0d8bp-1, -0x1.c73b23f5651bcp+0, 0x1.3195d5378b443p+2, -0x1.bd6451be65d02p+0, -0x1.abe3633c8ca4ap-1,
-0x1.ac615d2ace7a9p-1, -0x1.bd6451be65dp+0, 0x1.35d339eff7f5bp+2, -0x1.b8a1b55377a8p+0, -0x1.ab9db1ca214dcp-1,
-0x1.abe3633c8ca4ap-1, -0x1.b8a1b55377a8p+0, 0x1.381f85632e783p+2, -0x1.b5cf5e9643424p+0, -0x1.ab716cf0ec48ep-1,
-0x1.ab9db1ca214dap-1, -0x1.b5cf5e9643424p+0, 0x1.3991fd4a6bcfap+2, -0x1.b3f0a1ea0b17ep+0, -0x1.ab52ceb6f3f2ap-1,
-0x1.ab716cf0ec48ep-1, -0x1.b3f0a1ea0b17ep+0, 0x1.3a91119184b5ep+2, -0x1.b29bce180bf73p+0, -0x1.ab3c5dfbb87c8p-1,
-0x1.ab52ceb6f3f29p-1, -0x1.b29bce180bf73p+0, 0x1.3b4b862783e31p+2, -0x1.b19cb4dc8adadp+0, -0x1.ab2b36933f18cp-1,
-0x1.ab3c5dfbb87c8p-1, -0x1.b19cb4dc8adacp+0, 0x1.3bd9d02c84c06p+2, -0x1.b68f74bf26109p-1, 0,
-0x1.ab2b36933f18cp-1, -0x1.b68f74bf26109p-1, 0x1.a85ee1ffcc784p+2, 0, 0}, 10, 3);
/*
//Internal values for the Small Hamiltonian
std::vector<real> HSmall = {0x1.0426c96aab069p+0, 0x1.ded59bb2d5178p+0, 0x1p-51, -0x1.7p-49, 0x1.16cp-46,
0x1.ded59bb2d517bp+0, 0x1.55d0fc466a06dp+2, 0x1.9184bfa04e74ep+0, -0x1.cp-48, 0x1.2f9p-45,
0x0p+0, 0x1.9184bfa04e745p+0, 0x1.0c2e9d4f61f2fp+2, 0x1.fca8c90358575p-1, -0x1.d24p-48,
0x0p+0, 0x0p+0, 0x1.fca8c90358568p-1, 0x1.9ebfd12fa632p+2, 0x1.28638d909b166p+0,
0x0p+0, 0x0p+0, 0x0p+0, 0x1.28638d909b189p+0, 0x1.ea2b6ee3c0fcep+0};

//Internal values for Q Dagger
std::vector<real> QDag = {0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2, 0x1.43d136248490fp-2,
0x1.0fbbfac61672fp-1, -0x1.0107419ff28f8p-6, -0x1.f0941d57e1304p-3, -0x1.c4ec512b725b7p-3, -0x1.adc94797a50b3p-3, -0x1.9f58b15f3f51dp-3, -0x1.9573193949a4bp-3, -0x1.8e3c0c400880ap-3, 0x1.71e46a9daad3dp-4, 0x1.53a41530164a7p-1,
0x1.c225433b3398dp-8, -0x1.59dbbddb5928bp-2, -0x1.96a9ffb153486p-3, 0x1.3c676bd2640abp-2, 0x1.88911fe2eb707p-2, 0x1.613c0c90cd58cp-2, 0x1.5e1a6366fb9b2p-3, -0x1.d0d2b76cfc316p-2, -0x1.ca2d9a6cde297p-2, 0x1.c7cd610eab713p-3,
-0x1.2dc5d7bc76322p-2, -0x1.cffb69af0f358p-6, 0x1.2954196eb4af5p-10, 0x1.d60a82d12df5cp-3, -0x1.04df49758e38p-1, -0x1.0f4e4caa1638ep-2, 0x1.59c4e1985dddcp-1, -0x1.21802461169c9p-5, -0x1.443ee43aacb5dp-5, 0x1.10d1f64d50caep-2,
0x1.0269d9633fef1p-3, 0x1.c21bcf9e28034p-2, 0x1.cb6c3fbbc6ad9p-2, 0x1.8248a126305dbp-2, -0x1.f2f0ab70f7ffp-5, -0x1.428aa26ba2625p-2, -0x1.c74fa7c5e68e5p-3, -0x1.fc5340b545343p-2, -0x1.901ed6c7386b5p-3, -0x1.a0499570a38ep-4,
0x1.9239242ee51c2p-6, -0x1.57bd1e96bac59p-5, -0x1.14ca609f4623dp-2, 0x1.2ec9946cd0c03p-4, 0x1.1f48cf4797d4dp-3, -0x1.0db3bcba62a3fp-4, 0x1.79eea572eac44p-3, -0x1.05111bce60089p-1, 0x1.752fcc497ccf9p-1, -0x1.0e801f10e757ep-2};*/

void PrintArrayAsVec(std::string Name, la::band<real> & Array)
{
    std::cout << "std::vector<real> " << Name << " = {";
    int Sz = Array.NumElem();
    for (int i = 0; i < Sz-1; i++)
        std::cout << Array(i) << ", ";
    std::cout << Array(Sz-1) << "};" << std::endl;
}

std::vector<real> DivX2Compare = {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.4b91e11b56edbp-2, 0x1.483de460cbe8ep-3, 0x1.1e43d7edc9aa3p-5, 0x1.89360305d5e24p-9, 0x1.4e0effa870e26p-14, 0x1.89e46dbbb2165p-22, 0x1.4ebb6972f59b7p-35, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.483de460cbe8ep-3, 0x1.8393aca8ae411p-3, 0x1.9af44ea1a051ep-4, 0x1.7a0d23fce6d33p-6, 0x1.0eeaca711b803p-9, 0x1.dc77d6dcd714bp-15, 0x1.20f65549c5a91p-22, 0x1.f6b18e7b7a63p-36, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.1e43d7edc9aa3p-5, 0x1.9af44ea1a051ep-4, 0x1.fe7abc179c529p-4, 0x1.1a1bd36138786p-4, 0x1.0c9a70582612bp-6, 0x1.8c37e524befafp-10, 0x1.64ffd6e3e792dp-15, 0x1.ba0cb08530cc5p-23, 0x1.874b149ffa90cp-36, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.89360305d5e24p-9, 0x1.7a0d23fce6d33p-6, 0x1.1a1bd36138786p-4, 0x1.6a2b0902de681p-4, 0x1.9bb6319775d9ap-5, 0x1.919667ec07e1dp-7, 0x1.2e7086c55dbc7p-10, 0x1.157c2c5cd0d65p-15, 0x1.5d0b6f8be8304p-23, 0x1.3935479b410cfp-36, 0x0p+0, 0x0p+0, 0x1.4e0effa870e26p-14, 0x1.0eeaca711b803p-9, 0x1.0c9a70582612bp-6, 0x1.9bb6319775d9ap-5, 0x1.0e80230ae8c09p-4, 0x1.39d7e34ebdfb1p-5, 0x1.37abd7429b331p-7, 0x1.dcf21225efb52p-11, 0x1.bbc62f7db32b2p-16, 0x1.1a97fbdc9edc6p-23, 0x1.005d51f7e0844p-36, 0x0p+0, 0x1.89e46dbbb2165p-22, 0x1.dc77d6dcd714bp-15, 0x1.8c37e524befafp-10, 0x1.919667ec07e1dp-7, 0x1.39d7e34ebdfb1p-5, 0x1.a3a104c160f5bp-5, 0x1.ee7acd9403b3ep-6, 0x1.f1ecbf079a12fp-8, 0x1.81c259c2a7c4dp-11, 0x1.6af7dbab6680ep-16, 0x1.d2ef66a948438p-24, 0x0p+0, 0x1.4ebb6972f59b7p-35, 0x1.20f65549c5a91p-22, 0x1.64ffd6e3e792dp-15, 0x1.2e7086c55dbc7p-10, 0x1.37abd7429b331p-7, 0x1.ee7acd9403b3ep-6, 0x1.4f10dfe64f2d4p-5, 0x1.8faf0b39ed0cep-6, 0x1.96f32899fe82ep-8, 0x1.3e76cd146757ep-11, 0x1.2e6619ca469b7p-16, 0x0p+0, 0x0p+0, 0x1.f6b18e7b7a63p-36, 0x1.ba0cb08530cc5p-23, 0x1.157c2c5cd0d65p-15, 0x1.dcf21225efb52p-11, 0x1.f1ecbf079a12fp-8, 0x1.8faf0b39ed0cep-6, 0x1.11c4b2d02b5ap-5, 0x1.49cdbd358a973p-6, 0x1.52d95dcce5cd7p-8, 0x1.0b60da54aa975p-11, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.874b149ffa90cp-36, 0x1.5d0b6f8be8304p-23, 0x1.bbc62f7db32b2p-16, 0x1.81c259c2a7c4dp-11, 0x1.96f32899fe82ep-8, 0x1.49cdbd358a973p-6, 0x1.c7cffc81f09cep-6, 0x1.14cccb4a140d3p-6, 0x1.1e8a1fd7fa45p-8, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.3935479b410cfp-36, 0x1.1a97fbdc9edc6p-23, 0x1.6af7dbab6680ep-16, 0x1.3e76cd146757ep-11, 0x1.52d95dcce5cd7p-8, 0x1.14cccb4a140d3p-6, 0x1.81624ab952f02p-6, 0x1.d7464efb81c14p-7, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.005d51f7e0844p-36, 0x1.d2ef66a948438p-24, 0x1.2e6619ca469b7p-16, 0x1.0b60da54aa975p-11, 0x1.1e8a1fd7fa45p-8, 0x1.d7464efb81c14p-7, 0x1.4a20158d34538p-6, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0};
std::vector<real> HCompare = {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.7fe93e93e93e8p+1, 0x1.8d3d7e8292d35p-2, -0x1.5fb1706c5c1b2p+0, -0x1.e77a5b33ec226p-2, -0x1.30e5b9063b0eap-5, -0x1.0b27f850f96dp-11, -0x1.0cfeb60f94b02p-22, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.8d3d7e8292d35p-2, 0x1.7fe93e93e93e6p+1, 0x1.8d3d7e8292d3cp-2, -0x1.5fb1706c5c1afp+0, -0x1.e77a5b33ec223p-2, -0x1.30e5b9063b0e8p-5, -0x1.0b27f850f96cbp-11, -0x1.0cfeb60f94b05p-22, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.5fb1706c5c1b2p+0, 0x1.8d3d7e8292d3cp-2, 0x1.7fe93e93e93e8p+1, 0x1.8d3d7e8292d59p-2, -0x1.5fb1706c5c1b4p+0, -0x1.e77a5b33ec22ep-2, -0x1.30e5b9063b0ep-5, -0x1.0b27f850f96bep-11, -0x1.0cfeb60f94b03p-22, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.e77a5b33ec226p-2, -0x1.5fb1706c5c1afp+0, 0x1.8d3d7e8292d59p-2, 0x1.7fe93e93e93ebp+1, 0x1.8d3d7e8292d0ep-2, -0x1.5fb1706c5c1bp+0, -0x1.e77a5b33ec219p-2, -0x1.30e5b9063b0ep-5, -0x1.0b27f850f96cap-11, -0x1.0cfeb60f94b02p-22, 0x0p+0, 0x0p+0, -0x1.30e5b9063b0eap-5, -0x1.e77a5b33ec223p-2, -0x1.5fb1706c5c1b4p+0, 0x1.8d3d7e8292d0ep-2, 0x1.7fe93e93e93ebp+1, 0x1.8d3d7e8292d68p-2, -0x1.5fb1706c5c1b1p+0, -0x1.e77a5b33ec22bp-2, -0x1.30e5b9063b0e5p-5, -0x1.0b27f850f96ccp-11, -0x1.0cfeb60f94b34p-22, 0x0p+0, -0x1.0b27f850f96dp-11, -0x1.30e5b9063b0e8p-5, -0x1.e77a5b33ec22ep-2, -0x1.5fb1706c5c1bp+0, 0x1.8d3d7e8292d68p-2, 0x1.7fe93e93e93ecp+1, 0x1.8d3d7e8292d0ep-2, -0x1.5fb1706c5c1b1p+0, -0x1.e77a5b33ec224p-2, -0x1.30e5b9063b0f1p-5, -0x1.0b27f850f96dcp-11, 0x0p+0, -0x1.0cfeb60f94b02p-22, -0x1.0b27f850f96cbp-11, -0x1.30e5b9063b0ep-5, -0x1.e77a5b33ec219p-2, -0x1.5fb1706c5c1b1p+0, 0x1.8d3d7e8292d0ep-2, 0x1.7fe93e93e93eap+1, 0x1.8d3d7e8292d59p-2, -0x1.5fb1706c5c1b2p+0, -0x1.e77a5b33ec228p-2, -0x1.30e5b9063b0edp-5, 0x0p+0, 0x0p+0, -0x1.0cfeb60f94b05p-22, -0x1.0b27f850f96bep-11, -0x1.30e5b9063b0ep-5, -0x1.e77a5b33ec22bp-2, -0x1.5fb1706c5c1b1p+0, 0x1.8d3d7e8292d59p-2, 0x1.7fe93e93e93e7p+1, 0x1.8d3d7e8292d54p-2, -0x1.5fb1706c5c1b4p+0, -0x1.e77a5b33ec22ep-2, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.0cfeb60f94b03p-22, -0x1.0b27f850f96cap-11, -0x1.30e5b9063b0e5p-5, -0x1.e77a5b33ec224p-2, -0x1.5fb1706c5c1b2p+0, 0x1.8d3d7e8292d54p-2, 0x1.7fe93e93e93e7p+1, 0x1.8d3d7e8292d1ap-2, -0x1.5fb1706c5c1bp+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.0cfeb60f94b02p-22, -0x1.0b27f850f96ccp-11, -0x1.30e5b9063b0f1p-5, -0x1.e77a5b33ec228p-2, -0x1.5fb1706c5c1b4p+0, 0x1.8d3d7e8292d1ap-2, 0x1.7fe93e93e93ebp+1, 0x1.8d3d7e8292d51p-2, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.0cfeb60f94b34p-22, -0x1.0b27f850f96dcp-11, -0x1.30e5b9063b0edp-5, -0x1.e77a5b33ec22ep-2, -0x1.5fb1706c5c1bp+0, 0x1.8d3d7e8292d51p-2, 0x1.7fe93e93e93e4p+1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0};

std::vector<std::pair<real, real> > GaussQuadVal = {{-0x1.efb2b2ebf2106p-1, 0x1.4ce65f803eef8p-4}, {-0x1.ac0c44f0d0298p-1, 0x1.71f7a9b222beap-3}, {-0x1.3a0bd2077fd8cp-1, 0x1.0add87c827506p-2}, {-0x1.4c0916e48aa66p-2, 0x1.3fd7e9838d512p-2}, {0x0p+0, 0x1.522a43f65486ap-2}, {0x1.4c0916e48aa66p-2, 0x1.3fd7e9838d512p-2}, {0x1.3a0bd2077fd8cp-1, 0x1.0add87c827506p-2}, {0x1.ac0c44f0d0298p-1, 0x1.71f7a9b222beap-3}, {0x1.efb2b2ebf2106p-1, 0x1.4ce65f803eef8p-4}};

TEST(BSpline, SymmOverlapMatrix)    //Shortcut for symmetric matrices
{
    std::vector<real> Knots(18);
    real dx = 0.1;
    int k = 7;
    for (size_t i = 0; i < Knots.size(); i++)
        Knots[i] = i*dx;

    size_t No = Knots.size() - k;
//     size_t Nw = 2*k-1;

    quadrature::gauss<real, real> Gauss(GaussQuadVal);

    la::band<real> DivX2(No, k), H(No, k);
//     la::band<real> DivX2(No, Nw), H(No, Nw);
    spline::SymmOverlap(Gauss, Knots, k, DivX2, [](real x) {return (x ? 1.0 / (x*x) : 0.0);}, spline::BSpline<real>,  spline::BSpline<real>);
    spline::SymmOverlap(Gauss, Knots, k, H, [&Knots, k](int i, int j, real x) { return spline::DBSpline(k, i, x, Knots) * spline::DBSpline(k, j, x, Knots);});

    SCOPED_TRACE("Compare DivX2\n");
    Compare(DivX2, DivX2Compare);
    SCOPED_TRACE("Compare H\n");
    Compare(H, HCompare);
}

TEST(BSpline, OverlapMatrix)    //Shortcut for symmetric matrices
{
    std::vector<real> Knots(18);
    real dx = 0.1;
    int k = 7;
    for (size_t i = 0; i < Knots.size(); i++)
        Knots[i] = i*dx;

    size_t No = Knots.size() - k;
//     size_t Nw = 2*k-1;

    quadrature::gauss<real, real> Gauss(GaussQuadVal);

    la::band<real> DivX2(No, k), H(No, k);
    spline::Overlap(Gauss, Knots, k, DivX2, [](real x) {return (x ? 1.0 / (x*x) : 0.0);}, spline::BSpline<real>,  spline::BSpline<real>);
    spline::Overlap(Gauss, Knots, k, H, [&Knots, k](int i, int j, real x) { return spline::DBSpline(k, i, x, Knots) * spline::DBSpline(k, j, x, Knots);});

    SCOPED_TRACE("Compare DivX2\n");
    Compare(DivX2, DivX2Compare);
    SCOPED_TRACE("Compare H\n");
    Compare(H, HCompare);
}
#include <cmath>
TEST(Quadrature, Gaussian)
{
    quadrature::gauss<real, real> Gauss(9);
    real Test = Gauss.Quad(real(0), real(1), [](real x) -> real { return std::sin(x); }); //I'm sure there is a better way of picking the correct overload.
    quadrature::gauss<real, real> CompGauss(GaussQuadVal);
    real Test2 = CompGauss.Quad(real(0), real(1), [](real x) -> real { return std::sin(x); } );

    ASSERT_DOUBLE_EQ(Test, Test2) << "Gauss quad Differs: " << std::endl;

    real Test3 = 0x1.d6bafe095f2eap-2;
    ASSERT_DOUBLE_EQ(Test, Test3) << "Gauss quad final result Differs: " << std::endl;
}

/* **********************************************************************
 * 
 *                          Basic Linear Algebra Tests
 * 
 * **********************************************************************/
std::vector<real> RefVec = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0};
real RefDot = 4900;
//This test will see if the basic aspects of the la::slice type are working correctly.
TEST(LinearAlgebra, Slice)
{
    SCOPED_TRACE("Slice test\n");
    la::slice<real> A(RefVec);
    for (size_t i = 0; i < A.Size(); i++)
        ASSERT_DOUBLE_EQ(A[i], RefVec[i]) << "Slice initialisation-access failed.\n";
    la::slice<real> B(A);
    for (size_t i = 0; i < B.Size(); i++)
        ASSERT_DOUBLE_EQ(B[i], RefVec[i]) << "Slice copy constructor failed.\n";

    real Dot1 = Dot(A, A);
    real Dot2 = Dot(A, B);
    ASSERT_DOUBLE_EQ(Dot1, RefDot) << "Dot product failed.\n";
    ASSERT_DOUBLE_EQ(Dot2, RefDot) << "Dot product failed.\n";
}

//This test will see if the basic aspects of the la::vec type are working correctly.
TEST(LinearAlgebra, Vector)
{
    SCOPED_TRACE("Vector test\n");
    la::vec<real> A(RefVec.size());
    for (size_t i = 0; i < A.Size(); i++)
        A(i) = RefVec[i];
    la::vec<real> B(A);
    real Dot1 = Dot(A, A);
    real Dot2 = Dot(A, B);

    for (size_t i = 0; i < A.Size(); i++)
        ASSERT_DOUBLE_EQ(A(i), RefVec[i]) << "Slice initialisation-access failed.\n";
    for (size_t i = 0; i < B.Size(); i++)
        ASSERT_DOUBLE_EQ(B(i), RefVec[i]) << "Slice copy constructor failed.\n";
    ASSERT_DOUBLE_EQ(Dot1, RefDot) << "Dot product failed.\n";
    ASSERT_DOUBLE_EQ(Dot2, RefDot) << "Dot product failed.\n";
}


//TODO: Add unit tests for arrays
//TODO: Add a unit test making sure a diagonal array is the same as a k=1 banded array.

/*
 *
 * Krylov subspace test (Arnoldi method with a symmetric matrix)
 *
 */
TEST(Krylov, Arnoldi)
{
    SCOPED_TRACE("Krylov test\n");
    la::sqrarray<real> sqrH(1);
    sqrH.AddBlock(0, 0, &Ham);
    la::arnoldi<real> Kry(sqrH, 5);
    std::vector<real> Values = Kry.Eigenvalues();
    std::vector<real> Test = {0x1.da0fd16173623p-3, 0x1.9cec1088a4741p+0, 0x1.a322e7257174fp+1, 0x1.9ccce9d03f414p+2, 0x1.d7ea18c8e078ap+2};

//     io::Print(Values);
//     io::Print(Test);
    Compare(Test, Values);
}

int main(int argc, char **argv)
{
    std::cout << argc << " " << argv[0] << std::endl;
    try
    {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
    catch(ErrorCode i)
    {
        switch (i)
        {
            case UNKNOWN :
                std::cout << "UNKNOWN\n";
            break;
            case SIZE_MISMATCH :
                std::cout << "SIZE_MISMATCH\n";
            break;
            case DIM_MISMATCH :
                std::cout << "DIM_MISMATCH\n";
            break;
            case BLOCK_MISMATCH :
                std::cout << "BLOCK_MISMATCH\n";
            break;
            case OUT_OF_BOUNDS :
                std::cout << "OUT_OF_BOUNDS\n";
            break;
            case SELF_ASSIGNMENT :
                std::cout << "SELF_ASSIGNMENT\n";
            break;
        }
    }
//     return RUN_ALL_TESTS();
}
/*
void premain() __attribute__ ((constructor));
void premain()
{
        fputs("premain\n", stdout);
        char * a = "./newunittests";
        char ** b = &a;
        main(1, b);

}*/
