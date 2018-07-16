#include <vector>

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

nbody NBody;

typedef struct
{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} phase;

phase Phase;

vector<phase> Particles;
vector<phase> HaloParticles;
vector<phase> BulgeParticles;
vector<vector<phase> > Disk_Particles;
vector<vector<phase> > GasDisk_Particles;
