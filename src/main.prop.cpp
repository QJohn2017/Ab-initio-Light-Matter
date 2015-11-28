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
#include <iostream>
#include <map>
#include <libconfig.h++>
#include "laser/sine.h"
#include "laser/gauss.h"
#include "laser/trape.h"
#include "numeric/sequence.h"
#include "numeric/type.h"
#include "la/array.h"
#include "netcdf/get.h"

using namespace cathal;
using std::cout;
using std::endl;

void ConfigLaserSine(libconfig::Setting & Conf, laser::field<real, real> & Laser)
{
    real CEP, W0, E0, Tau, Shift;
    unsigned int Train, Cycles = 2;
    laser::carrier Shape = sin;

    Conf.lookupValue("CEP", CEP);
    Conf.lookupValue("PhotonEnergy", W0);
    Conf.lookupValue("MaxField", E0);
    Conf.lookupValue("Tau", Tau);
    Conf.lookupValue("Shift", Shift);
    Conf.lookupValue("NumTrains", Train);
    Conf.lookupValue("Cycles", Cycles);
//     Conf.lookup("Shape",&);

    Laser.AddPulse(new laser::sine(Train, Tau, Shift, E0, W0, CEP, Cycles, Shape));
}

void ConfigLaser(libconfig::Setting & Conf, laser::field<real, real> & Laser)
{
    libconfig::Setting & Sine = Conf.lookup("Sine"); //Exception check
    ConfigLaserSine(Sine, Laser);
}

void Config(libconfig::Config & Conf, laser::field<real, real> & Laser)
{
    libconfig::Setting & CLaser = Conf.lookup("Laser"); //Exception check
    ConfigLaser(CLaser, Laser);
}

int main(int argc, char * argv[])
{
    std::string CFile(argc > 1 ? argv[1] : "settings.laser.cfg");
/*****************************************************************
 *
 *          Load Data from NetCDF File.
 *
 * **************************************************************/
    std::string Dir = "in/";
    std::string File = Dir + "Basis.nc";
    std::string KName = "Knots";
    std::string EName = "Energy1d";
    std::string CName = "Coefficients";
    std::string DVName = "DipoleVelocity";
    std::string DLName = "DipoleMoment";

    std::string SName = "SplineOrder";
    std::string RName = "Radius";

//     la::fullblock<real> DipoleVelocity, DipoleMoment;

    std::vector<real> Knots = nc::GetVector<real>(File, KName);
    std::vector<real> Energy = nc::GetVector<real>(File, EName);
    std::cout << "Knots = " << Knots.size() << std::endl;

    la::fullblock<real> DipoleTest = nc::GetBlock<real>("in/Basis.0.1.nc", "DipoleMoment");

//     nc::Get(File, Energy, EName);
//     nc::Get(File, Knots, CName);
//     nc::Get(File, Knots, DVName);
//     nc::Get(File, Knots, DLName);
    unsigned int GaussN = 9;

    libconfig::Config InternalConf; //For settings you want to change less regularly.
    InternalConf.readFile(CFile.c_str());
    InternalConf.lookupValue("QuadOrder", GaussN);

    la::sqrarray<real> Hamiltonian(2);

    Hamiltonian.AddBlock(0, 1, &DipoleTest);
////////////////////////////////////////////////////////////////

//////TEST//////
    std::vector<size_t> In = {DipoleTest.Row(), DipoleTest.Column()};
    la::vec<real> Vec1(In), Vec2(In);
    Vec2.Set(1.0);
    Vec1 = Hamiltonian * Vec2;
///////////////
    libconfig::Config Conf;
    Conf.readFile(CFile.c_str());

    laser::field<real, real> Laser(GaussN);
    Config(Conf, Laser);

    int NumSteps = 1000;
    sequences::linear seq(Laser.End, NumSteps);
    real t = 0.0, Time = 0.0;
    while (!seq.End())
    {
        t = seq.Next();
        real E = Laser.E(t);
        real A = Laser.A(Time, t);
         cout << t << " " << E << " " << A << endl;

        Time = t;
    }
}

