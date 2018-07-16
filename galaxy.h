#ifndef _galaxy_h_
#define _galaxy_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>
#include <pthread.h>
#include <csignal>
#include <sstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>

#define PI M_PI

#define fourpi 12.56637061435917295385
#define twopi 6.28318530717958647693
#define oneover4pi 0.07957747154594766788
#define oneover2pi 0.15915494309189533577
#define oneoverpi 0.31830988618379067154
#define sqrt4pi 3.5449077018110320546
#define sqrt2pi 2.50662827463100050242
#define sqrtpi 1.7724538509055160273
#define oneoversqrt4pi 0.28209479177387814347
#define oneoversqrt2pi 0.39894228040143267794
#define oneoversqrtpi 0.56418958354775628695
#define sqrt2 1.4142135623730950488
#define oneoversqrt2 0.7071067811865475244
#define oneoversqrt8pi2 0.03582244801567226643

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

extern galaxy G;

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

extern astro_param Astro;

//Set up random number generator
//default is mt19937
extern const gsl_rng_type *T;
extern gsl_rng *rand_gen;
//gsl_rng_env_setup();
    
extern gsl_interp_accel *disk_mass_acc, *disk_mass_prime_acc;
extern gsl_interp_accel *bulge_mass_acc, *bulge_mass_prime_acc;
extern gsl_spline *disk_mass_spline;
extern gsl_spline *disk_mass_prime_spline;
extern gsl_spline *bulge_mass_spline;
extern gsl_spline *bulge_mass_prime_spline;
    
extern int n_psi, n_int, nr, nr_spline, n_iter, nrmx, nr_limit, lmax, l_max, max_iter, iter, n_simpson;
extern int halo_flag, disk_flag, gasdisk_flag, bulge_flag, smbh_flag, do_file_io, nbody_flag, chisq_flag, sersic_flag;
extern int nr_ac, contraction_flag, contraction_prescription;
extern int disk, gas_disk, gasdisk_params, disk_params, nondisk_params, astro_params, error_params;

extern double potcor, psi_crit, psi_d, psi_0, tolerance_factor, fraction;
extern double psi_upper_initial, psi_lower_initial, gamma_comp, log_rj2;
extern double dr, r_max, log_rmax, log_dr, delta_logr, baryon_frac;
extern double total_mass, halo_mass, disk_mass, total_disk_mass, bulge_mass;
extern double r_edge, halo_edge, disk_edge, bulge_edge;
extern double halo_const, bulge_const, fcut_halo, fcut_bulge;
extern double Q_total, X_total;
extern double sin_inclination, cos_inclination;

extern double Plcon[100];

//Use vectors instead of arrays because they are safer, just as fast (in principle) 
//and have better memory management, minimizing the chance of running out of memory
extern vector<double> Disk_Params, Disk_Const, Disk_Mass, Disk_Edge;
extern vector<double> GasDisk_Params, GasDisk_Const, GasDisk_Mass, GasDisk_Edge;
extern vector<double> Rho_Disk_Const, Rho_GasDisk_Const, Radius, S_1, S_2;
extern vector<double> Polytrope_Gas_Const;
extern vector<double> Table_E, DF_Sersic, DF_NFW, Dens_Psi_Bulge, Dens_Psi_Halo;
extern vector<double> H_Dens, H_Pot, H_FR, B_Dens, B_Pot, B_FR;
extern vector<double> D_Dens, D_Pot, D_FR, G_Dens, G_Pot, G_FR;
extern vector<double> Pot_Major_Tot, Pot_Minor, VC_Tot, Psi2, Pot_Major, Pot_Up;
extern vector<double> VC_Major, Vert_Freq, VC_Bulge, Vert_Freq_Bulge, Surf_Den;
extern vector<double> Omega, Omega2, A_K, A_K2;
extern vector<double> Am_Tab, R_Tab, R_Tab2, R_Tab2_Zero;
extern vector<double> Q_scale, X_scale, Q_avg, X_avg;
extern vector<double> Q_scale_gas, X_scale_gas, Q_avg_gas, X_avg_gas;

extern vector<double> Halo_Original_Radius, Halo_AC_Radius, Halo_AC_Dens; 
extern vector<double> Halo_AC_Dens_D, Halo_AC_Dens_DD, Halo_Mass_Radius;
extern vector<double> Bulge_Mass_Radius, Bulge_Mass_Prime;
extern vector<double> Disk_Mass_Radius, Disk_Mass_Prime;
extern vector<double> Raw_Halo_Dens, Raw_Halo_Dens_D, Raw_Halo_Dens_DD; 

extern vector<vector<double> > A_Dens;
extern vector<vector<double> > A_Pot;
extern vector<vector<double> > F_R;
extern vector<vector<double> > FR2;
extern vector<vector<double> > R_R;

extern vector<vector<double> > Halo_Dens;
extern vector<vector<double> > Halo_Pot;
extern vector<vector<double> > Halo_FR;

extern vector<vector<double> > Bulge_Dens;
extern vector<vector<double> > Bulge_Pot;
extern vector<vector<double> > Bulge_FR;

extern vector<vector<double> > Disk_Dens;
extern vector<vector<double> > Disk_Pot;
extern vector<vector<double> > Disk_FR;

extern vector<vector<double> > GasDisk_Dens;
extern vector<vector<double> > GasDisk_Pot;
extern vector<vector<double> > GasDisk_FR;

extern vector<vector<double> > Polytrope_Const;
extern vector<vector<double> > GasDensity_Const;
extern vector<vector<double> > GasDensity_Height;

extern vector<vector<double> > Rad_Spline, Rad_Spline_Gas;
extern vector<vector<double> > FD_Rat, D_Rat2, FSZ_Rat, SZ_Rat2;
extern vector<vector<double> > D_Rat, DZ2_Rat, FZ_Rat;
extern vector<vector<double> > FD_Rat_Gas, D_Rat2_Gas, FSZ_Rat_Gas, SZ_Rat2_Gas;
extern vector<vector<double> > D_Rat_Gas, DZ2_Rat_Gas, FZ_Rat_Gas;
extern vector<vector<double> > dirat;

extern vector<vector<vector<double> > > App_GasDisk_Pot, App_GasDisk_Force, App_GasDisk_Dens;
extern vector<vector<vector<double> > > App_GasDisk_Force_R, App_GasDisk_Force_Z;
extern vector<vector<vector<double> > > DPsi_DR, DPsi_DZ, D2Psi_DR2, D2Psi_DZ2;
extern vector<vector<vector<double> > > App_Disk_Pot, App_Disk_Force, App_Disk_Dens;
extern vector<vector<vector<double> > > App_Disk_Force_R, App_Disk_Force_Z;

extern vector<vector<double> > Legendre_Table;

//Function declarations

void   AllocateVectors(void);
void   AppDiskForce(double &s, double &z, double &fsad, double &fzad);
double AppDiskPot(double &s, double &z);
double AppDiskPot3(double &s, double &z);
double AppDiskDens(double &s, double &z);
double AppDiskDens3(double &s, double &z);
double AppGasDiskDens(double &s, double &z);
void   AppGasDiskForce(double &s, double &z, double &fsad, double &fzad);
double AppGasDiskPot(double &s, double &z);
double BulgeDens(double r, double z);
double BulgeDensPsi(double &energy);
double BulgeDF(double &energy);
void   BulgeForce(double &s, double &z, double &fs, double &fz);
void   BulgeMassRPrime(void);
void   BulgeMassWithRadius(void);
void   BulgePotential(void);
void   BulgePotentialEstimate(void);
void   CalculateGasDensConstants(void);
void   CalculateGasScaleHeights(void);
void   CalculatePolytropeConstants(void);
void   ClearVectors(void);
double coef(const int j);
void   ContractHalo(void);
void   DBH(void);
double Dens(double r, double z);
double DensrPsi(double &rad, double &psi);
double DF(double &energy);
double DiskDens(double r, double z, double psi);
double DiskDensEstimate(double s, double z);
double DiskDensf(double r, double z);
double DiskDensfI(double r, double z, int i);
double DiskDensI(double r, double z, double psi, int i);
double DiskDensity(double &r);
double DiskDensPsi(double &r, double &z, double &psi);
void   DiskDF(void);
double DiskDF3ez(double &ep, double &am, double &ez, int &i);
double DiskDF3ezGas(double &ep, double &am, double &ez, int &i);
double DiskDF3intez(double &ep, double &am, double &ez, int &i);
double DiskDF3intezGas(double &ep, double &am, double &ez, int &i);
double DiskDF5ez(double &vr, double &vt, double &vz, double &r, double &z, int &i);
double DiskDF5ezGas(double &vr, double &vt, double &vz, double &r, double &z, int &i);
double DiskDF5intez(double &vt, double &r, double z, int &i);
double DiskDF5intezGas(double &vt, double &r, double z, int &i);
double DiskForce(double &r);
void   DiskMassRPrime(void);
void   DiskMassWithRadius(void);
void   DiskPotentialEstimate(void);
void   DiskProfile(double &radius, double &z, int &i, 
                   double *Sigma_Profile, double *Rho_Profile);
void   DiskProfilePrime(double &radius, double &z, int &i, 
                        double *Sigma_Profile, double *Rho_Profile);
void   DiskProfile2Prime(double &radius, double &z, int &i, 
                         double *Sigma_Profile, double *Rho_Profile);
double dPolarDiskDens(double r, double cos_theta);
double dPolarGasDiskDens(double r, double cos_theta);
void   FindBrackets(double &psi, double &r_lower, double &r_upper);
double FnaMidDen(double &r, int &j);
double FnaMidDenGas(double &r, int &j);
void   Force(double &s, double &z, double &fs, double &fz);
double GasDiskDens(double r, double z, double psi);
double GasDiskDens2(double r, double z, double psi);
double GasDiskDensEstimate(double s, double z);
double GasDiskDensf(double r, double z);
double GasDiskDensf2(double r, double z);
double GasDiskDensfI(double r, double z, int i);
double GasDiskDensfI2(double r, double z, int i);
double GasDiskDensfI2(double r, double z, int i, double poly_const, 
                      double dens_const);
double GasDiskDensfI(double r, double z, int i, double poly_const);
double GasDiskDensI(double r, double z, double psi, int i, double poly_const);
double GasDiskDensI2(double r, double z, double psi, int i, 
                     double poly_const, double dens_const);
double GasDiskDensINoTrunc(double r, double z, double psi, int i, double poly_const);
double GasDiskDensINoTrunc2(double r, double z, double psi, int i, double dens_const, double psi_zero);
double GasDiskDensity(double &r);
double GasDiskDensPsi(double r, double z, double psi);
double GasDiskForce(double &r);
void   GasDiskPotentialEstimate(void);
void   GasDiskProfile(double &radius, double &z, int &i, 
                      double *Sigma_Profile, double *Rho_Profile);
void   GasDiskProfilePrime(double &radius, double &z, int &i,
                           double *Sigma_Profile, double *Rho_Profile);
void   GasDiskProfile2Prime(double &radius, double &z, int &i, 
                            double *Sigma_Profile, double *Rho_Profile);
double GasDiskSurfaceDensfI(double r, int j, double poly_const);
double GasDiskSurfaceDensfI(double r, int j, double poly_const, double zmax);
double GasDiskSurfaceDensfI2(double r, int j, double dens_const, double zmax);
double GasDiskSurfaceDensfI2(double r, int j, double dens_const);
void   GenBlackHole(double &mass);
void   GenDensPsiBulge(void);
void   GenDensPsiHalo(void);
double GetBulgeMass(double r);
double GetBulgeMassPrime(double r);
double GetDiskMassR(double r);
double GetDiskMassPrime(double r);
double GetGasHeight(double r, int j);
void   GetLegendreTable(void);
double GetNFWd2rhodpsi2(double &r);
void   GenNFWDistFunc(void);
void   GenSersicDistFunc(void);
void   GenTableE(void);
double Get_b_n(void);
double Get_rho_0(void);
void   GetAppDiskPot(void);
void   GetAppGasDiskPot(void);
double GetBulgeDens(double &psi);
double GetBulgePsi(double &r);
double GetDiskMass(int i);
void   GetDiskParameters(void);
double GetGasDiskMass(int i);
void   GetGasDiskParameters(void);
double GetDiskPsi(double &r, double &z);
void   GetEpicycleFrequencies(void);
void   GetFreqs(void);
double GetHaloConst(void);
double GetHaloDens(double &psi);
double GetHaloPsi(double &r);
double GetHaloTrunc(double &r);
double GetHaloTrunc2Prime(double &r);
double GetHaloTruncPrime(double &r);
void   GetNBodyRealisations(void);
double GetNFWDistFunc(double &energy);
void   GetOmegaKappa(double &r, double &freq_omega, double &freq_kappa);
void   GetParameters(void);
double GetPolyConst(double r, int j);
double GetDensConst(double r, int j);
void   GetRCirc(void);
void   GetRotationCurve(void);
double GetSersicd2rhodpsi2(double &r);
double GetSersicDistFunc(double &energy);
double GetTotalPsi(double &r);
double GetTrunc(double &r, int &i);
double GetTruncGas(double &r, int &i);
void   GetTruncPrime(double &r, int &i, double &truncfac, double &truncfacprime);
void   GetTruncPrimeGas(double &r, int &i, double &truncfac, double &truncfacprime);
double HaloDens(double r, double z);
double HaloDensity(double &r);
double HaloDensity2Prime(double &r);
double HaloDensityPrime(double &r);
void   HaloDensityProfile(double &r, double *profile);
double HaloDensPsi(double &energy);
double HaloDF(double energy);
double HaloForce(double &r);
void   HaloForce(double &s, double &z, double &fs, double &fz);
void   HaloPotential(void);
void   HaloPotentialEstimate(void);
double IterateGasConst(double r, int j);
double IterateGasDens(double r, int j);
double HaloProfileDens(double &r);
double HaloProfileDens2Prime(double &r);
double HaloProfileDensPrime(double &r);
double PolarBulgeDens(double r, double cos_theta, int l);
double PolarDens(double r, double cos_theta, int l);
double PolarHaloDens(double r, double cos_theta, int l);
double Pot(double s, double z);
double PotVector(double s, double z, vector<vector<double> >& Vec);
double RCirc(double &am);
double RawHaloProfile(double &r);
double RawHaloProfilePrime(double &r);
double RawHaloProfile2Prime(double &r);
void   ReadDensPsiBulge(void);
void   ReadDensPsiHalo(void);
void   ReadFreqs(void);
double RootBisection(double &psi, double &r_upper, double &r_lower, double tolerance);
double SersicDens(double &r);
double SersicDens2Prime(double &r);
double SersicDensPrime(double &r);
double SersicDensProfile(double &r, double *profile);
double SersicForce(double &r);
double SersicMass(double &r);
double SersicPotential(double &r);
double SigR2(double r, int &i);
double SigR2Gas(double r, int &i);
double SigR2Total(double r);
double SigZ2(double r, int &i);
double SigZ2Gas(double r, int &i);
double SigZ2Total(double r);
void   SplineD(vector<double> &x, vector<double> &y, int n, double yp1, 
               double ypn, vector<double> &y2);
void   SplintD(vector<double> &xa, vector<double> &ya, vector<double> &y2a, 
               int n, double x, double &y);
double TotalDens(double &r, double &z);
void   WriteBDat(void);
void   WriteDBHDat(void);
void   WriteHDat(void);

void   CopyGlobalsToDeviceChiSquare(void);
void   CopyGlobalsToDevice1(void);
void   CopyGlobalsToDevice2(void);
void   DiskPotentialEstimateCUDA(void);
void   GasDiskCUDA(void);
void   GetADensCUDA(const int l);
double GetSurfaceBrightnessCUDA(double error_factor, vector<double>& Radii, 
                                vector<double>& Data, vector<double>& Error, 
                                vector<double>& Ellipticity);
double GetCircVelocityCUDA(double error_factor, vector<double>& Radii, 
                           vector<double>& Data, vector<double>& Error);
double GetStarVelocityCUDA(double error_factor, vector<double>& Radii, 
                           vector<double>& Data, vector<double>& Error);
double GetStarDispersionCUDA(double error_factor, vector<double>& Radii, 
                             vector<double>& Data, vector<double>& Error);

void GenNFWDistFuncCUDA(void);
void GenSersicDistFuncCUDA(void);

#endif











