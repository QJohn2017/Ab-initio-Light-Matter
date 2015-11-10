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
typedef double pulsetype;
typedef double vecpot;

#include "laser/sine.h"
#include "laser/gauss.h"
#include "laser/trape.h"
#include "util/array.h"
#include "util/banded.h"
#include "util/io.h"
#include "numeric/sequence.h"
#include "numeric/type.h"
#include "numeric/krylov.h"
#include "numeric/splines.h"
#include <gtest/gtest.h>

using namespace cathal;
using std::cout;
using std::endl;

void PrintLoop(laser::field & Laser)
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

std::vector<std::array<real, 3>> Loop(laser::field & Laser)
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
    EXPECT_EQ(x.size(), y.size());   
    for (int i = 0; i < x.size(); i++)
    {
        EXPECT_EQ(x[i][0], y[i][0]);
        EXPECT_EQ(x[i][1], y[i][1]);
        EXPECT_EQ(x[i][2], y[i][2]);
    }
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
    laser::field Laser(GaussN);

    Laser.AddPulse(new laser::sine(train, 0.0, 0.0, cep, w0, e0, cycles, shape));
//     PrintLoop(Laser);
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
    unsigned int cycles = 10;
    laser::carrier shape = sin;
    unsigned int GaussN = 9;

    laser::field Laser(GaussN);

    Laser.AddPulse(new laser::gauss(train, tau, shift, cep, w0, e0, 10.0, 0.0, 5.0));
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

    laser::field Laser(GaussN);

    Laser.AddPulse(new laser::ctrape(1, 0.0, 2.0*shift, cep, w0, e0, 2.0, 10.0, shape));
//      PrintLoop(Laser);
    std::vector<std::array<real, 3>> Values = Loop(Laser);

    std::vector<std::array<real, 3>> Test {{0x0p+0, 0x0p+0, 0x0p+0},{0x1.3a142d88879c9p+3, 0x0p+0, 0x0p+0},{0x1.3a142d88879c9p+4, 0x0p+0, 0x0p+0},{0x1.d71e444ccb6aep+4, -0x1.c2e79a36b2739p-9, 0x1.598a7c84e6e6ap-3},{0x1.3a142d88879c9p+5, 0x1.46c377ddcdabdp-4, -0x1.92a05d90c1128p-3},{0x1.889938eaa983bp+5, -0x1.25f76ba8cafbfp-3, 0x1.023512a3f1166p-3},{0x1.d71e444ccb6aep+5, 0x1.7c5c8244b1613p-3, -0x1.66003a03aff64p-4},{0x1.12d1a7d776a9p+6, -0x1.99944c7104eb8p-3, -0x1.f2413614b0e18p-7},{0x1.3a142d88879c9p+6, 0x1.793a6855ee701p-3, 0x1.09261c8c052f1p-4},{0x1.6156b33998902p+6, -0x1.202bcabf2c337p-3, -0x1.3e1b0cde51f66p-3},{0x1.889938eaa983bp+6, 0x1.142a225dd0acdp-4, 0x1.3f12f85d7a258p-3},{0x1.afdbbe9bba774p+6, 0x0p+0, 0x1.a8b42f186f8p-13}};
     Compare(Values, Test);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
