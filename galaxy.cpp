#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace::std;

//Galaxy with all the input parameters plus b and p for the bulge
typedef struct
{
    double c_halo;
    double v_halo;
    double a_halo;
    double drtrunc_halo;
    double cusp;
    double halo_stream;
    double n_sersic;
    double ppp;
    double b_n;
    double v_bulge;
    double a_bulge;
    double bulge_stream;
	double rho_0;
    double bh_mass;
    vector<double> M_Disk;
    vector<double> R_Disk;
    vector<double> Out_Disk;
    vector<double> Z_Disk;
    vector<double> Dr_Trunc;
    vector<double> R_Kormendy;
    vector<double> Alpha;
    vector<double> Sigma_0;
    vector<double> R_Sigma;
    vector<double> Extra_Disk_Parameters;
    vector<double> Extra_DiskDF_Parameters;
    vector<double> M_GasDisk;
    vector<double> R_GasDisk;
    vector<double> Out_GasDisk;
    vector<double> Z_GasDisk;
    vector<double> Gamma;
    vector<double> Dr_Trunc_Gas;
    vector<double> R_Kormendy_Gas;
    vector<double> Alpha_Gas;
    vector<double> Sigma_0_Gas;
    vector<double> R_Sigma_Gas;
    vector<double> Extra_GasDisk_Parameters;
    vector<double> Extra_GasDiskDF_Parameters;
} galaxy;

galaxy G;

//Extra input parameters that do not affect the galaxy's intrinsic properties
typedef struct
{
    vector<double> ML_Disk;
    double ml_bulge;
    double inclination;
    double distance;
    double v_sys;
    double r_sys;
} astro_param;

astro_param Astro;

//Set up random number generator
//default is mt19937
const gsl_rng_type *T;
gsl_rng *rand_gen;
//gsl_rng_env_setup();
    
gsl_interp_accel *disk_mass_acc, *disk_mass_prime_acc;
gsl_interp_accel *bulge_mass_acc, *bulge_mass_prime_acc;
gsl_interp_accel *halo_dens_acc, *halo_dens_prime_acc, *halo_dens_2prime_acc;
gsl_spline *disk_mass_spline;
gsl_spline *disk_mass_prime_spline;
gsl_spline *bulge_mass_spline;
gsl_spline *bulge_mass_prime_spline;
gsl_spline *halo_dens_spline;
gsl_spline *halo_dens_prime_spline;
gsl_spline *halo_dens_2prime_spline;
    
//T = gsl_rng_mt19937;
//rand_gen = gsl_rng_alloc(T);

int n_psi, n_int, nr, nr_spline, n_iter, nrmx, nr_limit, lmax, l_max, max_iter, iter, n_simpson;
int halo_flag, disk_flag, gasdisk_flag, bulge_flag, smbh_flag, do_file_io, nbody_flag, chisq_flag, sersic_flag;
int nr_ac, contraction_flag, contraction_prescription;
int disk, gas_disk, disk_params, gasdisk_params, nondisk_params, astro_params, error_params;

double potcor, psi_crit, psi_d, psi_0, tolerance_factor, fraction;
double psi_upper_initial, psi_lower_initial, gamma_comp, log_rj2;
double dr, r_max, log_rmax, log_dr, delta_logr, baryon_frac;
double total_mass, halo_mass, disk_mass, total_disk_mass, bulge_mass;
double r_edge, halo_edge, disk_edge, bulge_edge;
double halo_const, bulge_const, fcut_halo, fcut_bulge;
double Q_total, X_total;
double sin_inclination, cos_inclination;

double Plcon[100];

vector<double> Disk_Params, Disk_Const, Disk_Mass, Disk_Edge;
vector<double> GasDisk_Params, GasDisk_Const, GasDisk_Mass, GasDisk_Edge;
vector<double> Rho_Disk_Const, Rho_GasDisk_Const, Radius, S_1, S_2;
vector<double> Table_E, DF_Sersic, DF_NFW, Dens_Psi_Bulge, Dens_Psi_Halo;
vector<double> H_Dens, H_Pot, H_FR, B_Dens, B_Pot, B_FR;
vector<double> D_Dens, D_Pot, D_FR, G_Dens, G_Pot, G_FR;
vector<double> Pot_Major_Tot, Pot_Minor, VC_Tot, Psi2, Pot_Major, Pot_Up;
vector<double> VC_Major, Vert_Freq, VC_Bulge, Vert_Freq_Bulge, Surf_Den;
vector<double> Omega, Omega2, A_K, A_K2;
vector<double> Am_Tab, R_Tab, R_Tab2, R_Tab2_Zero;
vector<double> Q_scale, X_scale, Q_avg, X_avg;
vector<double> Q_scale_gas, X_scale_gas, Q_avg_gas, X_avg_gas;

vector<double> Halo_Original_Radius, Halo_AC_Radius, Halo_AC_Dens; 
vector<double> Halo_AC_Dens_D, Halo_AC_Dens_DD, Halo_Mass_Radius;
vector<double> Bulge_Mass_Radius, Bulge_Mass_Prime;
vector<double> Disk_Mass_Radius, Disk_Mass_Prime;
vector<double> Raw_Halo_Dens, Raw_Halo_Dens_D, Raw_Halo_Dens_DD; 

vector<vector<double> > A_Dens;
vector<vector<double> > A_Pot;
vector<vector<double> > F_R;
vector<vector<double> > FR2;
vector<vector<double> > R_R;

vector<vector<double> > Halo_Dens;
vector<vector<double> > Halo_Pot;
vector<vector<double> > Halo_FR;

vector<vector<double> > Bulge_Dens;
vector<vector<double> > Bulge_Pot;
vector<vector<double> > Bulge_FR;

vector<vector<double> > Disk_Dens;
vector<vector<double> > Disk_Pot;
vector<vector<double> > Disk_FR;

vector<vector<double> > GasDisk_Dens;
vector<vector<double> > GasDisk_Pot;
vector<vector<double> > GasDisk_FR;

vector<vector<double> > Polytrope_Const;
vector<vector<double> > GasDensity_Const;
vector<vector<double> > GasDensity_Height;

vector<vector<double> > Rad_Spline, Rad_Spline_Gas;
vector<vector<double> > FD_Rat, D_Rat2, FSZ_Rat, SZ_Rat2;
vector<vector<double> > D_Rat, DZ2_Rat, FZ_Rat;
vector<vector<double> > FD_Rat_Gas, D_Rat2_Gas, FSZ_Rat_Gas, SZ_Rat2_Gas;
vector<vector<double> > D_Rat_Gas, DZ2_Rat_Gas, FZ_Rat_Gas;
vector<vector<double> > dirat;

vector<vector<vector<double> > > App_GasDisk_Pot, App_GasDisk_Force, App_GasDisk_Dens;
vector<vector<vector<double> > > App_GasDisk_Force_R, App_GasDisk_Force_Z;
vector<vector<vector<double> > > DPsi_DR, DPsi_DZ, D2Psi_DR2, D2Psi_DZ2;
vector<vector<vector<double> > > App_Disk_Pot, App_Disk_Force, App_Disk_Dens;
vector<vector<vector<double> > > App_Disk_Force_R, App_Disk_Force_Z;

vector<vector<double> > Legendre_Table;
