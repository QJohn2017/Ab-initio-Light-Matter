#include <iostream>
typedef double pulsetype;
typedef double vecpot;

#include "laser/sine.h"
#include "laser/gauss.h"
#include "laser/trape.h"
#include "numeric/sequence.h"
#include "numeric/type.h"

using namespace cathal;
using std::cout;
using std::endl;

int main(void)
{
    real cep = 0.0;
    real w0 = 1.0;
    real e0 = 0.2;
    real tau = 0.0;
    real shift = 1.0;
    unsigned int train = 1;
    unsigned int cycles = 2;
    laser::carrier shape = sin;
    unsigned int GaussN = 3;

    laser::field<real, real> Laser(GaussN);

    Laser.AddPulse(new laser::sine(train, tau, shift, e0, w0, cep, cycles, shape));

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

