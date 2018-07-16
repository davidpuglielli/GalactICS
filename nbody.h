#ifndef _nbody_h_
#define _nbody_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <pthread.h>
#include <csignal>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define PI M_PI

using namespace::std;

typedef struct
{
    vector<int> N_Disk;
    vector<int> Random_Seed;
    vector<double> Disk_Particle_Mass;
    vector<int> N_GasDisk;
    vector<int> Gas_Random_Seed;
    vector<double> GasDisk_Particle_Mass;
    int n_halo;
    int random_halo;
    double halo_particle_mass;
    int n_bulge;
    int random_bulge;
    double bulge_particle_mass;
} nbody;

extern nbody NBody;

typedef struct
{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} phase;

extern phase Phase;

extern vector<phase> Particles;
extern vector<phase> HaloParticles;
extern vector<phase> BulgeParticles;
extern vector<vector<phase> > Disk_Particles;
extern vector<vector<phase> > GasDisk_Particles;

void   GetDiskParticles(vector<phase> &Particles, int q, int random_seed);
void   GetGasDiskParticles(vector<phase> &Particles, int q, int random_seed);
void   GetBulgeParticles(vector<phase> &Particles, int random_seed);
void   GetHaloParticles(vector<phase> &Particles, int random_seed);
void   GetNBodyRealisations(void);
double InvertU(double &u);

#endif
