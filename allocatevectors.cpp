//A separate function to allocate the space required by all the global vectors
//used in the code. Putting this into a separate function allows us to call it
//once at the start of either GenerateGalaxy or MCMC. Otherwise the memory
//allocation of so many vectors (especially the two-dimensional ones) at each
//call of DBH() will slow things down. (I haven't actually tested this, but at
//the very least calling this function shouldn't hurt performance.)

//This code assumes that such quantities as nr and lmax will not change over
//the course of a full MCMC run. If that is not the case, AllocateVectors()
//will need to be called at the appropriate time in MCMC.cpp.

//Here is where we also modify nr if needed and initialise numerical quantities
//from in.gendenspsi, and clear the smaller vectors that were implemented with
//push_back() instead of being preallocated.

#include "galaxy.h"
#include "nbody.h"

void AllocateVectors(void)
{
    if (nr%2==0)
    {
        ++nr;
    }
    
    Radius.resize(nr);
    S_1.resize(nr);
    S_2.resize(nr);
    Table_E.resize(n_psi); 
    DF_Sersic.resize(n_psi); 
    DF_NFW.resize(n_psi); 
    Dens_Psi_Bulge.resize(n_psi); 
    Dens_Psi_Halo.resize(n_psi);
    H_Dens.resize(nr); 
    H_Pot.resize(nr); 
    H_FR.resize(nr); 
    B_Dens.resize(nr); 
    B_Pot.resize(nr); 
    B_FR.resize(nr); 
    D_Dens.resize(nr); 
    D_Pot.resize(nr); 
    D_FR.resize(nr);
    G_Dens.resize(nr); 
    G_Pot.resize(nr); 
    G_FR.resize(nr);
    Pot_Major_Tot.resize(nr); 
    Pot_Minor.resize(nr); 
    VC_Tot.resize(nr); 
    Psi2.resize(nr); 
    Pot_Major.resize(nr); 
    Pot_Up.resize(nr);
    VC_Major.resize(nr); 
    Vert_Freq.resize(nr); 
    VC_Bulge.resize(nr); 
    Vert_Freq_Bulge.resize(nr); 
    Surf_Den.resize(nr);
    Omega.resize(nr); 
    Omega2.resize(nr); 
    A_K.resize(nr); 
    A_K2.resize(nr); 
    Rad_Spline.resize(disk, vector<double>(nr_spline+1,0));
    Rad_Spline_Gas.resize(gas_disk, vector<double>(nr_spline+1,0));
    Am_Tab.resize(nr); 
    R_Tab.resize(nr); 
    R_Tab2.resize(nr); 
    R_Tab2_Zero.resize(nr);
    FD_Rat.resize(disk, vector<double>(nr_spline+1,0)); 
    D_Rat2.resize(disk, vector<double>(nr_spline+1,0)); 
    FSZ_Rat.resize(disk, vector<double>(nr_spline+1,0)); 
    SZ_Rat2.resize(disk, vector<double>(nr_spline+1,0));
    D_Rat.resize(disk, vector<double>(nr_spline+1,0)); 
    DZ2_Rat.resize(disk, vector<double>(nr_spline+1,0)); 
    FZ_Rat.resize(disk, vector<double>(nr_spline+1,0)); 
    FD_Rat_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 
    D_Rat2_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 
    FSZ_Rat_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 
    SZ_Rat2_Gas.resize(gas_disk, vector<double>(nr_spline+1,0));
    D_Rat_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 
    DZ2_Rat_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 
    FZ_Rat_Gas.resize(gas_disk, vector<double>(nr_spline+1,0)); 

    A_Dens.resize(l_max, vector<double>(nr,0));
    A_Pot.resize(l_max, vector<double>(nr,0));
    F_R.resize(l_max, vector<double>(nr,0));
    FR2.resize(l_max, vector<double>(nr,0));
    R_R.resize(l_max, vector<double>(nr,0));

    Halo_Dens.resize(l_max, vector<double>(nr,0));
    Halo_Pot.resize(l_max, vector<double>(nr,0));
    Halo_FR.resize(l_max, vector<double>(nr,0));

    Bulge_Dens.resize(l_max, vector<double>(nr,0));
    Bulge_Pot.resize(l_max, vector<double>(nr,0));
    Bulge_FR.resize(l_max, vector<double>(nr,0));

    Disk_Dens.resize(l_max, vector<double>(nr,0));
    Disk_Pot.resize(l_max, vector<double>(nr,0));
    Disk_FR.resize(l_max, vector<double>(nr,0));
    
    GasDisk_Dens.resize(l_max, vector<double>(nr,0));
    GasDisk_Pot.resize(l_max, vector<double>(nr,0));
    GasDisk_FR.resize(l_max, vector<double>(nr,0));
    
    Polytrope_Const.resize(gas_disk, vector<double>(nr,0));
    GasDensity_Const.resize(gas_disk, vector<double>(nr,0));
    GasDensity_Height.resize(gas_disk, vector<double>(nr,0.5));
    
    dirat.resize(5, vector<double>(nr_spline+1,0));
    
    if (contraction_flag)
    {
        if (nr_ac%2==0)
        {
            ++nr_ac;
        }
    
        Halo_Original_Radius.resize(nr_ac);
        Halo_AC_Radius.resize(nr_ac);
        Halo_AC_Dens.resize(nr_ac);
        Halo_AC_Dens_D.resize(nr_ac);
        Halo_AC_Dens_DD.resize(nr_ac);
        Halo_Mass_Radius.resize(nr_ac);
        Bulge_Mass_Radius.resize(nr_ac);
        Bulge_Mass_Prime.resize(nr_ac);
        Disk_Mass_Radius.resize(nr_ac);
        Disk_Mass_Prime.resize(nr_ac);
    }
    
//     App_GasDisk_Pot.resize(gas_disk);
//     App_GasDisk_Force.resize(gas_disk);
//     App_GasDisk_Dens.resize(gas_disk);
//     App_GasDisk_Force_R.resize(gas_disk);
//     App_GasDisk_Force_Z.resize(gas_disk);
//     DPsi_DR.resize(gas_disk);
//     D2Psi_DR2.resize(gas_disk);
//     DPsi_DZ.resize(gas_disk);
//     D2Psi_DZ2.resize(gas_disk);
//     
//     for (int i = 0; i < App_GasDisk_Pot.size(); ++i)
//     {
//         App_GasDisk_Pot[i].resize(nr, vector<double>(nr, 0));
//         App_GasDisk_Force[i].resize(nr, vector<double>(nr, 0));
//         App_GasDisk_Dens[i].resize(nr, vector<double>(nr, 0));
//         App_GasDisk_Force_R[i].resize(nr, vector<double>(nr, 0));
//         App_GasDisk_Force_Z[i].resize(nr, vector<double>(nr, 0));
//         DPsi_DR[i].resize(nr, vector<double>(nr, 0));
//         D2Psi_DR2[i].resize(nr, vector<double>(nr, 0));
//         DPsi_DZ[i].resize(nr, vector<double>(nr, 0));
//         D2Psi_DZ2[i].resize(nr, vector<double>(nr, 0));
//     }
    
//     App_Disk_Pot.resize(disk);
//     App_Disk_Force.resize(disk);
//     App_Disk_Dens.resize(disk);
//     App_Disk_Force_R.resize(disk);
//     App_Disk_Force_Z.resize(disk);
//     DPsi_DR.resize(disk);
//     D2Psi_DR2.resize(disk);
//     DPsi_DZ.resize(disk);
//     D2Psi_DZ2.resize(disk);
//     
//     for (int i = 0; i < App_Disk_Pot.size(); ++i)
//     {
//         App_Disk_Pot[i].resize(nr, vector<double>(nr, 0));
//         App_Disk_Force[i].resize(nr, vector<double>(nr, 0));
//         App_Disk_Dens[i].resize(nr, vector<double>(nr, 0));
//         App_Disk_Force_R[i].resize(nr, vector<double>(nr, 0));
//         App_Disk_Force_Z[i].resize(nr, vector<double>(nr, 0));
//         DPsi_DR[i].resize(nr, vector<double>(nr, 0));
//         D2Psi_DR2[i].resize(nr, vector<double>(nr, 0));
//         DPsi_DZ[i].resize(nr, vector<double>(nr, 0));
//         D2Psi_DZ2[i].resize(nr, vector<double>(nr, 0));
//     }
    
    Raw_Halo_Dens.resize(nr);
    Raw_Halo_Dens_D.resize(nr);
    Raw_Halo_Dens_DD.resize(nr);
    
    T = gsl_rng_mt19937;
    rand_gen = gsl_rng_alloc(T);
}

void ClearVectors(void)
{
    G.M_GasDisk.clear();
    G.R_GasDisk.clear();
    G.Out_GasDisk.clear();
    G.Z_GasDisk.clear();
    G.Dr_Trunc_Gas.clear();
    G.Sigma_0_Gas.clear();
    G.R_Sigma_Gas.clear();
    G.R_Kormendy_Gas.clear();
    G.Alpha_Gas.clear();
    G.Extra_GasDisk_Parameters.clear();
    
    G.M_Disk.clear();
    G.R_Disk.clear();
    G.Out_Disk.clear();
    G.Z_Disk.clear();
    G.Dr_Trunc.clear();
    G.Sigma_0.clear();
    G.R_Sigma.clear();
    G.R_Kormendy.clear();
    G.Alpha.clear();
    G.Extra_Disk_Parameters.clear();
    
    Disk_Const.clear();
    Rho_Disk_Const.clear();
    Disk_Edge.clear();
    Disk_Mass.clear();
    GasDisk_Const.clear();
    Rho_GasDisk_Const.clear();
    GasDisk_Edge.clear();
    GasDisk_Mass.clear();
    Q_scale.clear();
    X_scale.clear();
    Q_avg.clear();
    X_avg.clear();
    Q_scale_gas.clear();
    X_scale_gas.clear();
    Q_avg_gas.clear();
    X_avg_gas.clear();
        
    NBody.N_Disk.clear();
    NBody.Random_Seed.clear();
    NBody.Disk_Particle_Mass.clear();
    
    gsl_spline_free(disk_mass_spline);
    gsl_spline_free(disk_mass_prime_spline);
    gsl_spline_free(bulge_mass_spline);
    gsl_spline_free(bulge_mass_prime_spline);
    gsl_spline_free(halo_dens_spline);
    gsl_spline_free(halo_dens_prime_spline);
    gsl_spline_free(halo_dens_2prime_spline);
    gsl_interp_accel_free(disk_mass_acc);
    gsl_interp_accel_free(disk_mass_prime_acc);
    gsl_interp_accel_free(bulge_mass_acc);
    gsl_interp_accel_free(bulge_mass_prime_acc);
    gsl_interp_accel_free(halo_dens_acc);
    gsl_interp_accel_free(halo_dens_prime_acc);
    gsl_interp_accel_free(halo_dens_2prime_acc);
}
