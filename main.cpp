#include "lattice.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
#include <fstream>

using namespace std;

int main()
{
    cout << "hello world" << endl;
    int nx = 250;
    int ny = 100;
    int nz = 40;
    vector<LBM::Lattice> lattices(nx*ny*nz, LBM::Lattice());
    LBM::Lattice::calRelaventParameters();
    LBM::Lattice::initializeLatticeProperties(lattices,nx,ny,nz);
    LBM::Lattice::readLatticeType4file(lattices);
    cout << "Initialization done..." << endl;
    cout << "Solid fraction = " << LBM::Lattice::checkSolidFraction(lattices) << endl;

    for (const auto& lattice : lattices) {
        for (int i = 0; i < LBM::vSetLength; ++i) {
            if (lattice.neighLattices[i] == nullptr) {
                cout << "neighbour lattice is not correctly set..." << endl;
            }
        }
    }


    clock_t timerTotal = clock();
    clock_t timer = clock();
    int infoInterval2screen = 100;
    int infoInterval2file = 10;
    int LongInfoInterval2file = 1000;

    LBM::Lattice::writeLatticeType2binaryfile(lattices);
    LBM::Lattice::writeBasicInfo2file(lattices);
    
    FILE* fileMomentumEnergy = fopen("momentumEnergy.dat", "w");

    for (int i = 0; i <= 5000; ++i) {
        LBM::Lattice::runCoreFunc(lattices);

        if (fmod(i,infoInterval2screen) == 0) {

            timer = clock() - timer;
            printf("Current step = %d, Time used per %d steps = %lf, Total time used = %lf\n", i, infoInterval2screen,
                (static_cast<double>(timer) / CLOCKS_PER_SEC),
                ((static_cast<double>(clock()) - timerTotal) / CLOCKS_PER_SEC));
            timer = clock();
        }

        if (fmod(i,infoInterval2file) == 0) {
            cout << "Mass conservation check: " << LBM::Lattice::checkMassConservation(lattices) << endl;
            cout << "Maximum velocity check: " << LBM::Lattice::checkMaximumVelocity(lattices) << endl;
            vector<double> me = LBM::Lattice::getMomentumEnergy(lattices);
            fprintf(fileMomentumEnergy, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", me[0], me[1], me[2], me[3], me[4], me[5], me[6]);
            fflush(fileMomentumEnergy);
        }

        if (fmod(i,LongInfoInterval2file) == 0) {
            LBM::Lattice::write2binaryfile(lattices, i);
        }
    }
    
    fclose(fileMomentumEnergy);

    return 0;
}
