//Generate the N-body realisations
//Input file in.nbody is of the form
//
// <halo_particles>   <halo_random_seed>
// <bulge_particles>  <bulge_random_seed>
// <disk_particles#1> <disk_random_seed#1>
// <disk_particles#2> <disk_random_seed#2>
// etc.
//
//for as many disks as found in in.diskpars.

#include "galaxy.h"
#include "chisq.h"
#include "nbody.h"

void GetNBodyRealisations(void)
{
    cout << "Generating N-Body realisations..." << endl;
    
    nbody_flag = 1;
    
    double time = 0;
    
    ifstream nbodyfile("in.nbody", ios::in);
    
    if (!nbodyfile.is_open())
    {
        cerr << "Can't find or can't open in.nbody. Please see README for more "
             << "info. Exiting..." << endl;
        exit(1);
    }
   
    //Get halo and bulge particle numbers 
    nbodyfile >> NBody.n_halo >> NBody.random_halo;
    NBody.halo_particle_mass = halo_mass/NBody.n_halo;
    
    cout << "Halo particles: " << NBody.n_halo << ", using seed = " 
         << NBody.random_halo << endl;
    cout << "Halo particle mass = " << NBody.halo_particle_mass << endl;
    
    nbodyfile >> NBody.n_bulge >> NBody.random_bulge;
    NBody.bulge_particle_mass = bulge_mass/NBody.n_bulge;
    
    cout << "Bulge particles: " << NBody.n_bulge << ", using seed = " 
         << NBody.random_bulge << endl;
    cout << "Bulge particle mass = " << NBody.bulge_particle_mass << endl;
    
    //Get particle numbers for the stellar disk(s)
    for (int i = 0; i < disk; ++i)
    {
        int nn;
        nbodyfile >> nn;
        
        if (nbodyfile.good())
        {
            NBody.N_Disk.push_back(nn);
        }
        else
        {
            cerr << "   in.nbody is not the correct form, or you have too "
                 << "few nbody\n   parameters for the disk. Exiting..." << endl;
            exit(1);
        }
        
        nbodyfile >> nn;
        
        if (nbodyfile.good())
        {
            NBody.Random_Seed.push_back(nn);
        }
        else
        {
            cerr << "   in.nbody is not the correct form, or you have too "
                 << "few nbody\n   parameters for the disk. Exiting..." << endl;
            exit(1);
        }
        
        NBody.Disk_Particle_Mass.push_back(Disk_Mass[i]/NBody.N_Disk[i]);
        
        cout << "Disk particles: " << NBody.N_Disk[i] << ", using seed = "
             << NBody.Random_Seed[i] << endl;
        cout << "Disk particle mass = " << NBody.Disk_Particle_Mass[i] << endl;
    }
    
    //Get particle numbers for the gas disk(s)
    for (int i = 0; i < gas_disk; ++i)
    {
        int nn;
        nbodyfile >> nn;
        
        if (nbodyfile.good())
        {
            NBody.N_GasDisk.push_back(nn);
        }
        else
        {
            cerr << "   in.nbody is not the correct form, or you have too "
                 << "few nbody\n   parameters for the gas disk. Exiting..." << endl;
            exit(1);
        }
        
        nbodyfile >> nn;
        
        if (nbodyfile.good())
        {
            NBody.Gas_Random_Seed.push_back(nn);
        }
        else
        {
            cerr << "   in.nbody is not the correct form, or you have too "
                 << "few nbody\n   parameters for the disk. Exiting..." << endl;
            exit(1);
        }
        
        NBody.GasDisk_Particle_Mass.push_back(GasDisk_Mass[i]/NBody.N_GasDisk[i]);
        
        cout << "Disk particles: " << NBody.N_GasDisk[i] << ", using seed = "
             << NBody.Gas_Random_Seed[i] << endl;
        cout << "Disk particle mass = " << NBody.GasDisk_Particle_Mass[i] << endl;
    }
    
    cout << "Read in nbody data" << endl;
    cout << endl;
    
    //Get the halo particles and write to halo file
    
    if (halo_flag)
    {
        HaloParticles.resize(NBody.n_halo);
    
        GetHaloParticles(HaloParticles, NBody.random_halo);
    
        ofstream halofile("halo", ios::out | ios::binary);
        
        //FILE *fhalo;
        //fhalo = fopen("halo1","w");
    
        halofile.write((char *)&NBody.n_halo, sizeof(int));
        halofile.write((char *)&time, sizeof(double));
        for (int i = 0; i < NBody.n_halo; ++i)
        {
            halofile.write((char *)&HaloParticles[i], sizeof(phase));
        }
        
        //fwrite(&NBody.n_halo,sizeof(int),1,fhalo);
        //fwrite(&time,sizeof(double),1,fhalo);
        //for(int i=0;i<NBody.n_halo;++i)
        //{
        //    fwrite(&HaloParticles[i].mass,sizeof(phase),1,fhalo);
        //}
        
        //ofstream halo2("halo2");
        //for (int i = 0; i < NBody.n_halo; ++i)
        //{
        //    halo2 << HaloParticles[i].mass << " " << HaloParticles[i].x
        //          << " " << HaloParticles[i].y << " " << HaloParticles[i].z << " "
        //          << HaloParticles[i].vx << " " << HaloParticles[i].vy << " "
        //          << HaloParticles[i].vz << endl;
        //}
    
        halofile.close();
        
        cout << "Halo particles written to file 'halo'." << endl;
    }
    
    //Get the bulge particles and write to bulge file
    
    if (bulge_flag)
    {
        BulgeParticles.resize(NBody.n_bulge);
    
        GetBulgeParticles(BulgeParticles, NBody.random_bulge);
    
        ofstream bulgefile("bulge", ios::out | ios::binary);
     
        bulgefile.write((char *)&NBody.n_bulge, sizeof(int));
        bulgefile.write((char *)&time, sizeof(double));
        for (int i = 0; i < NBody.n_bulge; ++i)
        {
            bulgefile.write((char *)&BulgeParticles[i], sizeof(phase));
        }
    
        bulgefile.close();
        
        cout << "Bulge particles written to file 'bulge'." << endl;
    }
    
    //Now get the disk particles sequentially for each component
    
    if (disk_flag)
    {
        Disk_Particles.resize(disk);

        for (int i = 0; i < disk; ++i)
        {
            Disk_Particles[i].resize(NBody.N_Disk[i]);

            GetDiskParticles(Disk_Particles[i], i, NBody.Random_Seed[i]);

            stringstream disk_number;

            disk_number << i+1;
            string filename = "disk" + disk_number.str();

            ofstream diskfile(filename.c_str(), ios::out | ios::binary);

            diskfile.write((char *)&NBody.N_Disk[i], sizeof(int));
            diskfile.write((char *)&time, sizeof(double));
            for (int j = 0; j < NBody.N_Disk[i]; ++j)
            {
                diskfile.write((char *)&Disk_Particles[i][j], sizeof(phase));
            }

            diskfile.close();
            
            cout << "Disk " << i+1 << " particles written to file '" 
                 << filename << "'" << endl;
        }
    }
    
    //Now get the gas disk particles sequentially for each component
    
    if (gasdisk_flag)
    {
        GasDisk_Particles.resize(gas_disk);

        for (int i = 0; i < gas_disk; ++i)
        {
            GasDisk_Particles[i].resize(NBody.N_GasDisk[i]);

            GetGasDiskParticles(GasDisk_Particles[i], i, NBody.Gas_Random_Seed[i]);

            stringstream disk_number;

            disk_number << i+1;
            string filename = "gasdisk" + disk_number.str();

            ofstream gasdiskfile(filename.c_str(), ios::out | ios::binary);

            gasdiskfile.write((char *)&NBody.N_GasDisk[i], sizeof(int));
            gasdiskfile.write((char *)&time, sizeof(double));
            for (int j = 0; j < NBody.N_GasDisk[i]; ++j)
            {
                gasdiskfile.write((char *)&GasDisk_Particles[i][j], sizeof(phase));
            }

            gasdiskfile.close();
            
            cout << "Gas Disk " << i+1 << " particles written to file '" 
                 << filename << "'" << endl;
        }
    }
}

void GetHaloParticles(vector<phase> &Particles, int random_seed)
{
    //double rtest = 0.0001097863074, ztest = -8.254986097e-05, psitest = 23.87755302;
    //cout << "halodens " << HaloDens(rtest,ztest) << " " << HaloDensPsi(psitest) << endl;
    
    //Set up random number generator
    //default is mt19937
    //const gsl_rng_type *T;
    //gsl_rng *rand_gen;
    //gsl_rng_env_setup();
    
    //T = gsl_rng_mt19937;
    //rand_gen = gsl_rng_alloc(T);
    gsl_rng_set(rand_gen, random_seed);
    
    //JD: `Find maximum of rho*r^2
    //DP: In *two* directions
    
    double rho_max = 0, rho_min;
    double u1_max = halo_edge, v1_max = 0.5*PI;
    double broadcon = 1;
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i], rho_temp, z = 0;
        //double rad = i*0.01*dr, rho_temp, z = 0;
        
        rho_temp = HaloDens(rad, z);
        rho_temp *= rad*rad;
        
        if (rho_temp > rho_max)
        {
            rho_max = rho_temp;
        }
        
        //cout << i << " " << rho_temp << endl;
    }
    
    cout << "rhomax " << rho_max << endl;
    
    for (int i = 1; i < nr; ++i)
    {
        double rad = 0, rho_temp, z = Radius[i];
        
        rho_temp = HaloDens(rad, z);
        rho_temp *= z*z;
        
        if (rho_temp > rho_max)
        {
            rho_max = rho_temp;
        }
    }
    
    cout << "rhomax " << rho_max << endl;
    
    rho_max *= 1.5;
    rho_min = 1e-10*rho_max;
    
    cout << "Calculating halo positions and velocities..." << endl;
    
    //Note: Calls to HaloDens produce different results from the original code. The
    //original code seems to give the right answer before this but the wrong answer
    //(much smaller for small R) when called in the following loop. Calling halodens
    //with radius and z arguments gives different results from calling halodenspsi
    //with psi argument using psi=pot(radius,z).
    
    for (int i = 0; i < NBody.n_halo; )
    {
        double u1 = u1_max*gsl_rng_uniform(rand_gen);
        double v1 = 2*v1_max*(gsl_rng_uniform(rand_gen)-0.5);
        
        double rad = u1;
        double z = rad*tan(v1);//cout<<"u  "<<v1<<" "<<rad<< " " << z << " " << ceil(sqrt(rad*rad+z*z)/dr) << endl;
//         if (z > halo_edge)
//         {
//             z = halo_edge;
//         }
        
        double rho_test = HaloDens(rad, z)*(rad*rad+z*z);
        
        if(rho_test < rho_min)
        {
            continue;
        }
        
        double rho_ran = (rho_max - rho_min)*gsl_rng_uniform(rand_gen);
        
        if (rho_ran > rho_test)
        {
            continue;
        }
        
        double phi = 2*PI*gsl_rng_uniform(rand_gen);
        double x = rad*cos(phi);
        double y = rad*sin(phi);
        double v_r, v_phi, v_z;
        
        double psi = Pot(rad, z);
        
        if(psi < psi_crit)
            continue;
        
        double v_max2 = 2*(psi-psi_crit);
        double v_max = sqrt(v_max2);
        
        int count = 0;
        double f0 = 0, frand = 1; //JD: dummy starters
        double f_max = HaloDF(psi)-fcut_halo;
        
        if(f_max>100)
            cout << "fmax  " << i << " " << f_max << " " << psi << endl;
        
        while(frand > f0)
        {
            double v2 = 1.1*v_max2, E;
            
            while(v2 > v_max2)
            {
                v_r   = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                v_phi = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                v_z   = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                
                v2 = v_r*v_r + v_phi*v_phi + v_z*v_z;
                E = psi-0.5*v2;
                //cout << E << endl;
            }
            
            f0 = HaloDF(E)-fcut_halo;
            frand = f_max*gsl_rng_uniform(rand_gen);
            ++count;
        }
        
        //cout << i << endl;
        
        if (gsl_rng_uniform(rand_gen) < G.halo_stream)
        {
            v_phi = fabs(v_phi);
        }
        else
        {
            v_phi = -fabs(v_phi);
        }
        
        double cos_phi = x/rad, sin_phi = y/rad;
        double v_x = v_r*cos_phi - v_phi*sin_phi;
        double v_y = v_r*sin_phi + v_phi*cos_phi;
        
        Particles[i].mass = NBody.halo_particle_mass;
        Particles[i].x = x;
        Particles[i].y = y;
        Particles[i].z = z;
        Particles[i].vx = v_x;
        Particles[i].vy = v_y;
        Particles[i].vz = v_z;
        
        ++i;
        
        if (i%1000==0)
        {
            cout << "." << flush;
        }
    }
    
    gsl_rng_free(rand_gen);
    cout << endl;
}
    
void GetBulgeParticles(vector<phase> &Particles, int random_seed)
{
    //Set up random number generator
    //default is mt19937
    //const gsl_rng_type *T;
    //gsl_rng *rand_gen;
    //gsl_rng_env_setup();
    
    T = gsl_rng_mt19937;
    rand_gen = gsl_rng_alloc(T);
    gsl_rng_set(rand_gen, random_seed);
    
    //JD: `Find maximum of rho*r^2
    //DP: In *two* directions
    
    double rho_max = 0, rho_min;
    double u1_max = bulge_edge, v1_max = 0.5*PI;
    double broadcon = 1;
    
    for (int i = 1; i < nr; ++i)
    {
        double rad = Radius[i], rho_temp, z = 0;
        
        rho_temp = BulgeDens(rad, z);
        rho_temp *= rad*rad;
        
        if (rho_temp > rho_max)
        {
            rho_max = rho_temp;
        }
    }
    
    for (int i = 1; i < nr; ++i)
    {
        double rad = 0, rho_temp, z = Radius[i];
        
        rho_temp = BulgeDens(rad, z);
        rho_temp *= z*z;
        
        if (rho_temp > rho_max)
        {
            rho_max = rho_temp;
        }
    }
    
    cout << "rhomax " << rho_max << endl;
    
    rho_max *= 1.5;
    rho_min = 1e-10*rho_max;
    
    cout << "Calculating bulge positions and velocities..." << endl;
    
    for (int i = 0; i < NBody.n_bulge; )
    {
        double u1 = u1_max*gsl_rng_uniform(rand_gen);
        double v1 = 2*v1_max*(gsl_rng_uniform(rand_gen)-0.5);
        
        double rad = u1;
        double z = rad*tan(v1);
        double rho_test = BulgeDens(rad, z)*(rad*rad+z*z);
        
        if(rho_test < rho_min)
            continue;
        
        double rho_ran = (rho_max - rho_min)*gsl_rng_uniform(rand_gen);
        
        if (rho_ran > rho_test)
            continue;
        
        double phi = 2*PI*gsl_rng_uniform(rand_gen);
        double x = rad*cos(phi);
        double y = rad*sin(phi);
        double v_r, v_phi, v_z;
        
        double psi = Pot(rad, z);
        
        if(psi < psi_crit)
            continue;
        
        double v_max2 = 2*(psi-psi_crit);
        double v_max = sqrt(v_max2);
        
        int count = 0;
        double f0 = 0, frand = 1; //JD: dummy starters
        double f_max = BulgeDF(psi)-fcut_bulge;
        
        while(frand > f0)
        {
            double v2 = 1.1*v_max2, E;
            
            while(v2 > v_max2)
            {
                v_r   = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                v_phi = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                v_z   = 2*v_max*(gsl_rng_uniform(rand_gen)-0.5);
                
                v2 = v_r*v_r + v_phi*v_phi + v_z*v_z;
                E = psi-0.5*v2;
            }
            
            f0 = BulgeDF(E)-fcut_bulge;
            frand = f_max*gsl_rng_uniform(rand_gen);
            ++count;
        }
        
        if (gsl_rng_uniform(rand_gen) < G.bulge_stream)
        {
            v_phi = fabs(v_phi);
        }
        else
        {
            v_phi = -fabs(v_phi);
        }
        
        double cos_phi = x/rad, sin_phi = y/rad;
        double v_x = v_r*cos_phi - v_phi*sin_phi;
        double v_y = v_r*sin_phi + v_phi*cos_phi;
        
        Particles[i].mass = NBody.bulge_particle_mass;
        Particles[i].x = x;
        Particles[i].y = y;
        Particles[i].z = z;
        Particles[i].vx = v_x;
        Particles[i].vy = v_y;
        Particles[i].vz = v_z;
        
        ++i;
        if (i%1000==0)
        {
            cout << "." << flush;
        }
    }
    
    gsl_rng_free(rand_gen);
    cout << endl;
}
    
void GetDiskParticles(vector<phase> &Particles, int q, int random_seed)
{
    //Set up random number generator
    //default is mt19937
    //const gsl_rng_type *T;
    //gsl_rng *rand_gen;
    //gsl_rng_env_setup();
    
    T = gsl_rng_mt19937;
    rand_gen = gsl_rng_alloc(T);
    gsl_rng_set(rand_gen, random_seed);
    
    double rdisk = 1.2*G.R_Disk[q], zdisk = 1.2*G.Z_Disk[q];
    double rho_max = 0, rho_min = 0;
    double broadcon = 1;
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i], z = 0;
        double rho_guess = exp(-rad/G.R_Disk[q]);
        double rho_test = DiskDensfI(rad, z, q)/rho_guess;
        
        if (rho_test > rho_max)
        {
            rho_max = rho_test;
        }
    }
    
    rho_max *= 1.2;
    
    cout << "Calculating disk positions and velocities... " << rho_max << endl;
    
    int density_trials = 0, velocity_trials = 0;
    
    for (int i=0; i < NBody.N_Disk[q]; )
    {
        double f_max = -1;
        double u1, v1;
        double rad, z, v_r, v_phi, v_z;
        double v_phi_max, vphimax_old;
        double vsig_r, vsig_p, vsig_z;
        double x, y, phi, v_x, v_y;
        double omega, kappa;
        
        while (f_max < 0)
        {
            rad = 2*Disk_Edge[q];
            
            while (rad > Disk_Edge[q])
            {
                u1 = -gsl_rng_uniform(rand_gen);
                v1 = 2*(gsl_rng_uniform(rand_gen)-0.5);
                rad = rdisk*InvertU(u1);
                z = zdisk*atanh(v1);
            }
        
            double zcon = cosh(z/zdisk);
            double rho_guess = exp(-rad/rdisk)/(zcon*zcon);
            double rho_test = DiskDensfI(rad, z, q)/rho_guess;

            ++density_trials;

            if (rho_test < rho_min)
                continue;

            double rho_ran = (rho_max-rho_min)*gsl_rng_uniform(rand_gen);

            if (rho_ran > rho_test)
                continue;

            //cout << "Still in after outer while loop   " << endl;
            
            phi = 2*PI*gsl_rng_uniform(rand_gen);
            x = rad*cos(phi);
            y = rad*sin(phi);

            GetOmegaKappa(rad, omega, kappa);

            v_phi_max = omega*rad;
            vsig_r = sqrt(SigR2(rad, q));
            vsig_p = kappa*vsig_r/(2*omega);
            vsig_z = sqrt(SigZ2(rad, q));
            vphimax_old = v_phi_max;

            FindMax(rad, z, vsig_p, v_phi_max, f_max, q);
            f_max *= 1.1;
        }
        
        
        v_phi_max = vphimax_old;
        
        double f0 = 0, frand = 1, gr, gp, gz, g2;
        
        while(frand > f0)
        {
            g2 = 999;
            
            while (g2 > 1.0)
            {
                gr = 8*(gsl_rng_uniform(rand_gen)-0.5);
                gp = 16*(gsl_rng_uniform(rand_gen)-0.5);
                gz = 8*(gsl_rng_uniform(rand_gen)-0.5);
                g2 = gr*gr/16+gp*gp/64+gz*gz/16;
            }
            
            v_r   = broadcon*vsig_r*gr;
            v_phi = v_phi_max+broadcon*vsig_p*gp;
            v_z   = broadcon*vsig_z*gz;
            
            f0 = DiskDF5ez(v_r, v_phi, v_z, rad, z, q);
            frand = f_max*gsl_rng_uniform(rand_gen);
            
//             if (f0 > f_max)
//             {
//                 cout << "There's some kind of error here. "
//                      << rad << " " << z << endl;
//             }
            
            ++velocity_trials;
        }
        
        v_phi_max = v_phi-broadcon*vsig_p*gp;
        v_phi = v_phi_max+broadcon*vsig_p*gp;
        
        double cos_phi = x/rad, sin_phi = y/rad;
        v_x = v_r*cos_phi - v_phi*sin_phi;
        v_y = v_r*sin_phi + v_phi*cos_phi;
        
        double angular_momentum = rad*v_phi;
        double r_circ = RCirc(angular_momentum);
        
        GetOmegaKappa(r_circ, omega, kappa);
        //double vsig_zz = sqrt(SigR2(r_circ, q))*kappa/2/omega;
        
        Particles[i].mass = NBody.Disk_Particle_Mass[q];
        Particles[i].x = x;
        Particles[i].y = y;
        Particles[i].z = z;
        Particles[i].vx = v_x;
        Particles[i].vy = v_y;
        Particles[i].vz = v_z;
        
        //cout<<i<<endl;
        ++i;
        
        if (i%1000==0)
        {
            cout << "." << flush;
        }
    }
    
    cout << endl;
    cout << "Number of density trials: " << density_trials << endl;
    cout << "Number of velocity trials: " << velocity_trials << endl;
    
    gsl_rng_free(rand_gen);
}
        
void GetGasDiskParticles(vector<phase> &Particles, int q, int random_seed)
{
    //Set up random number generator
    //default is mt19937
    //const gsl_rng_type *T;
    //gsl_rng *rand_gen;
    //gsl_rng_env_setup();
    
    T = gsl_rng_mt19937;
    rand_gen = gsl_rng_alloc(T);
    gsl_rng_set(rand_gen, random_seed);
    
    double rdisk = 1.2*G.R_GasDisk[q], zdisk = 1.2*G.Z_GasDisk[q];
    double rho_max = 0, rho_min = 0;
    double broadcon = 1;
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i], z = 0;
        double rho_guess = exp(-rad/G.R_GasDisk[q]);
        double rho_test = GasDiskDensfI(rad, z, q)/rho_guess;
        
        if (rho_test > rho_max)
        {
            rho_max = rho_test;
        }
    }
    
    rho_max *= 1.2;
    
    cout << "Calculating gas disk positions and velocities... " << rho_max << endl;
    
    int density_trials = 0, velocity_trials = 0;
    
    for (int i=0; i < NBody.N_GasDisk[q]; )
    {
        double f_max = -1;
        double u1, v1;
        double rad, z, v_r, v_phi, v_z;
        double v_phi_max, vphimax_old;
        double vsig_r, vsig_p, vsig_z;
        double x, y, phi, v_x, v_y;
        double omega, kappa;
        
        while (f_max < 0)
        {
            rad = 2*GasDisk_Edge[q];
            
            while (rad > GasDisk_Edge[q])
            {
                u1 = -gsl_rng_uniform(rand_gen);
                v1 = 2*(gsl_rng_uniform(rand_gen)-0.5);
                rad = rdisk*InvertU(u1);
                z = zdisk*atanh(v1);
            }
        
            double zcon = cosh(z/zdisk);
            double rho_guess = exp(-rad/rdisk)/(zcon*zcon);
            double rho_test = GasDiskDensfI2(rad, z, q)/rho_guess;
            
            ++density_trials;

            if (rho_test < rho_min)
                continue;

            double rho_ran = (rho_max-rho_min)*gsl_rng_uniform(rand_gen);

            if (rho_ran > rho_test)
                continue;

            phi = 2*PI*gsl_rng_uniform(rand_gen);
            x = rad*cos(phi);
            y = rad*sin(phi);

            GetOmegaKappa(rad, omega, kappa);

            v_phi_max = omega*rad;
            vsig_r = sqrt(SigR2Gas(rad, q));
            vsig_p = kappa*vsig_r/(2*omega);
            vsig_z = sqrt(SigZ2Gas(rad, q));
            vphimax_old = v_phi_max;

            FindMaxGas(rad, z, vsig_p, v_phi_max, f_max, q);
            f_max *= 1.1;
        }
        
        v_phi_max = vphimax_old;
        
        double f0 = 0, frand = 1, gr, gp, gz, g2;
        
        while(frand > f0)
        {
            g2 = 999;
            
            while (g2 > 1.0)
            {
                gr = 8*(gsl_rng_uniform(rand_gen)-0.5);
                gp = 16*(gsl_rng_uniform(rand_gen)-0.5);
                gz = 8*(gsl_rng_uniform(rand_gen)-0.5);
                g2 = gr*gr/16+gp*gp/64+gz*gz/16;
            }
            
            v_r   = broadcon*vsig_r*gr;
            v_phi = v_phi_max+broadcon*vsig_p*gp;
            v_z   = broadcon*vsig_z*gz;
            
            f0 = DiskDF5ezGas(v_r, v_phi, v_z, rad, z, q);
            frand = f_max*gsl_rng_uniform(rand_gen);
            
//             if (f0 > f_max)
//             {
//                 cout << "There's some kind of error here. "
//                      << rad << " " << z << endl;
//             }
            
            ++velocity_trials;
        }
        
        v_phi_max = v_phi-broadcon*vsig_p*gp;
        v_phi = v_phi_max+broadcon*vsig_p*gp;
        
        double cos_phi = x/rad, sin_phi = y/rad;
        v_x = v_r*cos_phi - v_phi*sin_phi;
        v_y = v_r*sin_phi + v_phi*cos_phi;
        
        double angular_momentum = rad*v_phi;
        double r_circ = RCirc(angular_momentum);
        
        GetOmegaKappa(r_circ, omega, kappa);
        //double vsig_zz = sqrt(SigR2(r_circ, q))*kappa/2/omega;
        
        Particles[i].mass = NBody.GasDisk_Particle_Mass[q];
        Particles[i].x = x;
        Particles[i].y = y;
        Particles[i].z = z;
        Particles[i].vx = v_x;
        Particles[i].vy = v_y;
        Particles[i].vz = v_z;
        
        //cout<<i<<endl;
        ++i;
        
        if (i%1000==0)
        {
            cout << "." << flush;
        }
    }
    
    cout << endl;
    cout << "Number of density trials: " << density_trials << endl;
    cout << "Number of velocity trials: " << velocity_trials << endl;
    
    gsl_rng_free(rand_gen);
}
        
double InvertU(double &u)
{
    //JD: invert the u function to find R
    
    double rg=1, drr=1, eps=1e-8, r_new;
    int i = 0;
    
    while (i<20 && drr>eps)
    {
        r_new = rg-(-(1+rg)*exp(-rg)-u)/(rg*exp(-rg));
        drr = fabs(r_new-rg);
        rg = r_new;
        ++i;
    }
    
    if (i==20)
    {
        cerr << "Warning: R did not converge";
    }
    
    return rg;
}
