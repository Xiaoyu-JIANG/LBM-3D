#pragma once

#include <iostream>
#include <vector>
#include <omp.h>
#include <random>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>

#define D3Q19

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

namespace LBM
{
    enum class TYPE {
        INIT = 0, FLUID = 1, SOLID = 2
    };

    constexpr int NUM_THREADS = 8;

    #ifdef D3Q27
        constexpr int vSetLength = 27;
        // D3Q27 model p.87
        // lattice weights
        constexpr double weight[vSetLength] = {8./27.,
                                            2./27., 2./27., 2./27., 2./27., 2./27., 2./27.,
                                            1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54.,
                                            1./216., 1./216.,1./216.,1./216.,1./216.,1./216.,1./216.,1./216.}; // lattice weights
        // velocity set direction
        //                                     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
        constexpr int    dirX[vSetLength]   = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
        constexpr int    dirY[vSetLength]   = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1};
        constexpr int    dirZ[vSetLength]   = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1};
        // reverse directions
        constexpr int    revdir[vSetLength] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25};
        // velocity set c_i
        double cx[vSetLength]               = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
        double cy[vSetLength]               = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1};
        double cz[vSetLength]               = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1};
    #endif

    #ifdef D3Q19
        constexpr int vSetLength = 19;
        // D3Q19 model p.87
        // lattice weights
        constexpr double weight[vSetLength] = {1./3.,
                                               1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
                                               1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.}; // lattice weights
        // velocity set direction
        //                                     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
        constexpr int    dirX[vSetLength]   = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
        constexpr int    dirY[vSetLength]   = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
        constexpr int    dirZ[vSetLength]   = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
        // reverse directions
        constexpr int    revdir[vSetLength] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17};
        // velocity set c_i
        double cx[vSetLength]               = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
        double cy[vSetLength]               = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
        double cz[vSetLength]               = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
    #endif

    #ifdef D3Q15
        constexpr int vSetLength = 15;
        // D3Q15 model p.87
        // lattice weights
        constexpr double weight[vSetLength] = {2./9.,
                                            1./9., 1./9., 1./9., 1./9., 1./9., 1./9.,
                                            1./72., 1./72.,1./72.,1./72.,1./72.,1./72.,1./72.,1./72.}; // lattice weights
        // velocity set direction
        //                                     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14
        constexpr int    dirX[vSetLength]   = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
        constexpr int    dirY[vSetLength]   = {0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1};
        constexpr int    dirZ[vSetLength]   = {0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1};
        // reverse directions
        constexpr int    revdir[vSetLength] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13};
        // velocity set c_i
        double cx[vSetLength]               = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
        double cy[vSetLength]               = {0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1};
        double cz[vSetLength]               = {0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1};
    #endif

    

    double tau = 0.933012701892219;
    double dt = 0.0;
    double omega = 0.0;
    double cs2 = 0.0;
    double cs2inv = 0.0;
    double cs4inv = 0.0;
    double dynVisco = 1.0016e-3;
    double ds = 0.0001; // Lattice size
    double vol = ds*ds*ds;
    double density = 1000.0; // Fluid reference density

    int nx = 0;
    int ny = 0;
    int nz = 0;

    

    class Lattice
    {
    public:
        Lattice() {}

        void printNeighbours() {
            std::cout << neighLattices[0] << std::endl;
        }

        static void runCoreFunc(std::vector<Lattice>& lattices)
        {

            #pragma omp parallel
            {
                #pragma omp for
                for (int idx = 0; idx < lattices.size(); ++idx) lattices[idx].swapPopulation();

                #pragma omp for
                for (int idx = 0; idx < lattices.size(); ++idx) lattices[idx].streaming();

                #pragma omp for
                for (int idx = 0; idx < lattices.size(); ++idx) lattices[idx].collision();

            }


            // #pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
            // for (int idx = 0; idx < lattices.size(); ++idx) {
            //     lattices[idx].swapPopulation();
            // }
            // #pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
            // for (int idx = 0; idx < lattices.size(); ++idx) {
            //     lattices[idx].streaming();
            // }
            // #pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
            // for (int idx = 0; idx < lattices.size(); ++idx) {
            //     lattices[idx].collision();
            // }
            return;
        }

        static void initializeLatticeProperties(std::vector<Lattice>& lattices, const int nx, const int ny, const int nz)
        {

            omp_set_num_threads(NUM_THREADS);

            LBM::nx = nx;
            LBM::ny = ny;
            LBM::nz = nz;

            std::default_random_engine generator(0);
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            for (int k = 0; k < nz; ++k) {
                double zpos = k * ds;
                for (int j = 0; j < ny; ++j) {
                    double ypos = j * ds;
                    for (int i = 0; i < nx; ++i) {
                        double xpos = i * ds;
                        int idx = k * nx * ny + j * nx + i;
                        auto& lattice = lattices[idx];
                        for (int m = 0; m < vSetLength; ++m) {
                            int ii = (i + dirX[m] + nx) % nx;
                            int jj = (j + dirY[m] + ny) % ny;
                            int kk = (k + dirZ[m] + nz) % nz;
                            int idx_nb = kk * nx * ny + jj * nx + ii;
                            lattice.neighLattices[m] = &lattices[idx_nb];
                        }
                        lattice.type = TYPE::FLUID;
                        lattice.dens = density;
                        lattice.id = idx;
                        lattice.vx = 0.0;
                        lattice.vy = 0.0;
                        lattice.vz = 0.0;
                        lattice.fx = 1.0;
                        lattice.fy = 0.0;
                        lattice.fz = 0.0;
                        lattice.x = xpos;
                        lattice.y = ypos;
                        lattice.z = zpos;
                        lattice.ix = i;
                        lattice.iy = j;
                        lattice.iz = k;
                        lattice.initializeEquilibrium();
                    }
                }
            }
            return;
        }

        static void calRelaventParameters()
        {
            const double A = 1.0718;
            tau = (density * ds * ds * (1.0 - A / 2.0)) / (dynVisco * 3.0 * A * A);
            dt = A * tau;
            std::cout << "calculated dt = " << dt << std::endl;
            omega = dt / tau;
            for (int i = 0; i < vSetLength; ++i) {
                cx[i] *= (ds/dt);
                cy[i] *= (ds/dt);
                cz[i] *= (ds/dt);
            }
            cs2 = ds * ds / (3.0 * dt * dt);
            cs2inv = 1.0 / cs2;
            cs4inv = cs2inv * cs2inv;
            return;
        }

        static double checkMassConservation(std::vector<Lattice>& lattices)
        {
            double mass = 0.0;
            #pragma omp parallel for num_threads(NUM_THREADS) schedule(static) reduction(+:mass)
            for (int idx = 0; idx < lattices.size(); ++idx) {
                mass += lattices[idx].dens * vol;
            }
            return mass;
        }

        static double checkSolidFraction(std::vector<Lattice>& lattices)
        {
            double ns = 0;
            for (const auto& lattice : lattices) {
                if (lattice.type == TYPE::SOLID) ns++;
            }
            return ns/lattices.size();
        }

        static double checkMaximumVelocity(std::vector<Lattice>& lattices)
        {
            double velmax = 0.0;
            for (const auto& lattice : lattices) {
                double vel = sqrt(SQ(lattice.vx)+SQ(lattice.vy)+SQ(lattice.vz));
                if (vel > velmax) velmax = vel;
            }
            return velmax;
        }

        static std::vector<double> getMomentumEnergy(std::vector<Lattice>& lattices)
        {
            std::vector<double> momentum(7,0.0);
            for (const auto& lattice : lattices) {
                if (lattice.type == TYPE::SOLID) continue;
                double mass = lattice.dens * vol;
                momentum[0] += mass;
                momentum[1] += mass * lattice.vx;
                momentum[2] += mass * lattice.vy;
                momentum[3] += mass * lattice.vz;
                momentum[4] += mass * lattice.vx * lattice.vx;
                momentum[5] += mass * lattice.vy * lattice.vy;
                momentum[6] += mass * lattice.vz * lattice.vz;
            }
            return momentum;
        }

        static void write2binaryfile(std::vector<Lattice>& lattices, const int timeStep)
        {
            std::stringstream filename;
            filename << "massVel_" << std::setw(10) << std::setfill('0') << timeStep << ".bin";
            FILE* file = fopen(filename.str().c_str(), "wb");
            for (const auto& lattice : lattices) {
                double mass = lattice.dens * vol;
                fwrite(&mass, sizeof(double), 1, file);
                fwrite(&lattice.vx, sizeof(double), 1, file);
                fwrite(&lattice.vy, sizeof(double), 1, file);
                fwrite(&lattice.vz, sizeof(double), 1, file);
            }
            fflush(file);
            fclose(file);
            return;
        }

        static void writeLatticeType2binaryfile(std::vector<Lattice>& lattices)
        {
            FILE* file = fopen("latticeType.bin", "wb");
            int isFluid = 1;
            int isSolid = 0;
            for (const auto& lattice : lattices) {
                if (lattice.type == TYPE::FLUID) fwrite(&isFluid, sizeof(int), 1, file);
                if (lattice.type == TYPE::SOLID) fwrite(&isSolid, sizeof(int), 1, file);
            }
            fflush(file);
            fclose(file);
            return;
        }

        static void writeBasicInfo2file(std::vector<Lattice>& lattices)
        {
            std::ofstream file("basicInfo.dat");
            file << dt << "\t\t\ttime interval\n";
            file << density << "\t\t\tfluid reference density\n";
            file << dynVisco << "\t\t\tfluid dynamic viscosity\n";
            file << tau << "\t\t\trelaxtion time\n";
            file << ds << "\t\t\tlattice size\n";
            file << checkSolidFraction(lattices) << "\t\t\tsolid fraction\n";
            file << LBM::nx << "\t" << LBM::ny << "\t" << LBM::nz << "\t\t\tgeometrical size\n";
            file.close();
        }

        static void readLatticeType4file(std::vector<Lattice>& lattices)
        {
            FILE* file = fopen("initSpace.dat", "r");
            int dummy;
            for (auto& lattice : lattices) {
                fscanf(file, "%d", &dummy);
                if (dummy == 1) lattice.type = TYPE::FLUID;
                else if (dummy == 0) lattice.type = TYPE::SOLID;
            }
            return;
        }

        // member variables.
        int id = -1;
        TYPE type = TYPE::INIT;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        int ix = 0;
        int iy = 0;
        int iz = 0;

        double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        double dens = 0.0;
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;

        #ifdef D3Q27
            Lattice* neighLattices[vSetLength] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
            double pop_cur[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            double pop_pre[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            double pop_eq[vSetLength]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            double forcing[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        #endif

        #ifdef D3Q19
            Lattice* neighLattices[vSetLength] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
            double pop_cur[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double pop_pre[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double pop_eq[vSetLength]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double forcing[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        #endif

        #ifdef D3Q15
            Lattice* neighLattices[vSetLength] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
            double pop_cur[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double pop_pre[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double pop_eq[vSetLength]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            double forcing[vSetLength] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        #endif

    private:
        void collision(void)
        {
            // Skip if the current lattice is a solid one
            if (type == TYPE::INIT || type == TYPE::SOLID) return;
            // Update density and velocity
            updateDensVelo();
            // Calculate the equilibrium state
            calEquilibrium();
            // Collision step
            collide();
            return;
        }

        void streaming(void)
        {
            // Skip solid lattice
            if (type == TYPE::INIT || type == TYPE::SOLID) return;
            for (int i = 0; i < vSetLength; ++i) {
                auto neib_cell = neighLattices[revdir[i]];
                if (neib_cell->type == TYPE::SOLID) {
                    // pop_cur[i] = pop_pre[revdir[i]];
                    // LatticeBoltzmannMethod.pdf p.180
                    double CiDotVel = cx[i] * neib_cell->vx + cy[i] * neib_cell->vy + cz[i] * neib_cell->vz;
                    pop_cur[i] = pop_pre[revdir[i]] + 2.0 * cs2inv * weight[revdir[i]] * dens * CiDotVel;
                }
                else if (neib_cell->type == TYPE::FLUID) {
                    pop_cur[i] = neib_cell->pop_pre[i];
                }
            }
            return;
        }

        void swapPopulation(void)
        {
            for (int i = 0; i < vSetLength; ++i) pop_pre[i] = pop_cur[i];
            return;
        }

        void initializeEquilibrium(void)
        {
            calEquilibrium();
            for (int i = 0; i < vSetLength; ++i) {
                pop_cur[i] = pop_eq[i];
                pop_pre[i] = pop_eq[i];
            }
            return;
        }

        // Collides populations with body force
        void collide(void)
        {
            calBodyForce();
            for (int i = 0; i < vSetLength; ++i) {
                pop_cur[i] = (1.0 - omega) * pop_cur[i] + omega * pop_eq[i] + forcing[i] * dt;
            }
            return;
        }

        // Calculate source term of body force
        void calBodyForce(void)
        {
            for (int i = 0; i < vSetLength; ++i) {
                double CiDotVel = (cx[i] * vx) + (cy[i] * vy) + (cz[i] * vz);
                forcing[i] = (1.0 - omega/2.0) * weight[i] * (fx * (cs2inv * (cx[i] - vx) + cs4inv * CiDotVel * cx[i])
                                                            + fy * (cs2inv * (cy[i] - vy) + cs4inv * CiDotVel * cy[i])
                                                            + fz * (cs2inv * (cz[i] - vz) + cs4inv * CiDotVel * cz[i]));
            }
            return;
        }

        // Calculate equilibrium
        void calEquilibrium()
        {
            double velsq = SQ(vx) + SQ(vy) + SQ(vz);
            for (int i = 0; i < vSetLength; ++i) {
                double CiDotVel = (cx[i] * vx) + (cy[i] * vy) + (cz[i] * vz);
                //pop_eq[i] = weight[i] * (dens + 3.0 * CiDotVel + 4.5 * SQ(CiDotVel) - 1.5 * velsq);
                // From "A new partial-bounceback lattice-Boltzmann method for fluid flow through heterogeneous media"
                // !! density term is outside the parentheses
                // https://www.sciencedirect.com/science/article/pii/S009830040800201X?casa_token=V3l6htc9CcUAAAAA:OAmwGW4vlzrTedfW5MLDzwmTnkWqgc4D8YLC_nKPZJ6kbYUHQOMYH1dzujq70FI3vadPNpbehoMG
                pop_eq[i] = weight[i] * dens * (1.0 + cs2inv * CiDotVel + 0.5 * cs4inv * SQ(CiDotVel) - 0.5 * cs2inv * velsq);
                // pop_eq[i] = weight[i] * (dens + 3.0 * CiDotVel + 4.5 * SQ(CiDotVel) - 1.5 * velsq);
            }
            return;
        }

        // Update density and fluid velocity
        void updateDensVelo(void)
        {
            dens = 0.0;
            vx = 0.0;
            vy = 0.0;
            vz = 0.0;
            for (int i = 0; i < vSetLength; ++i) {
                dens += pop_cur[i];
                vx += cx[i] * pop_cur[i];
                vy += cy[i] * pop_cur[i];
                vz += cz[i] * pop_cur[i];
            }
            vx /= dens;
            vy /= dens;
            vz /= dens;
            vx += fx / (2.0 * dens) * dt;
            vy += fy / (2.0 * dens) * dt;
            vz += fz / (2.0 * dens) * dt;
            return;
        }




        void setNeighbourSolid(const int numRadius, const int iter, const int ix0, const int iy0, const int iz0)
        {
            if (iter > 0) {
                for (auto& neig : neighLattices) {
                    if (sqrt(double((ix0-neig->ix)*(ix0-neig->ix)+
                                    (iy0-neig->iy)*(iy0-neig->iy)+
                                    (iz0-neig->iz)*(iz0-neig->iz))) < numRadius)
                    {
                        neig->type = TYPE::SOLID;
                    }
                    neig->setNeighbourSolid(numRadius, iter-1, ix0, iy0, iz0);
                }

            }
            else {
                return;
            }
        }

    };
}
