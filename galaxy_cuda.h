#ifndef _galaxycuda_h_
#define _galaxycuda_h_

#include <thrust/host_vector.h>
#include <thrust/for_each.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/random.h>
#include <thrust/iterator/counting_iterator.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "galaxy.h"
#include "chisq.h"

#define DISK_LIMIT 10

//put the global variables here? Update them on the host?
__constant__ int n_psi_d, n_int_d, nr_d, nr_spline_d, n_iter_d, nrmx_d, lmax_d, l_max_d; 
__constant__ int max_iter_d, iter_d, n_simpson_d;
__constant__ int halo_flag_d, disk_flag_d, gasdisk_flag_d, bulge_flag_d, smbh_flag_d;
__constant__ int do_file_io_d, nbody_flag_d, chisq_flag_d, sersic_flag_d;
__constant__ int disk_d, disk_params_d, nondisk_params_d;
__constant__ int gas_disk_d, gasdisk_params_d;
__constant__ int astro_params_d, error_params_d;
__constant__ int do_chisq_file_io_d, n_los_d, n_vel_d;
__constant__ int nr_ac_d, contraction_flag_d, contraction_prescription_d;

__constant__ double potcor_d, psi_crit_d, psi_d_d, psi_0_d, tolerance_factor_d, fraction_d;
__constant__ double psi_upper_initial_d, psi_lower_initial_d, gamma_comp_d, log_rj2_d;
__constant__ double dr_d, r_max_d, log_rmax_d, log_dr_d, delta_logr_d, a00_d;
__constant__ double total_mass_d, halo_mass_d, disk_mass_d, total_disk_mass_d, bulge_mass_d;
__constant__ double r_edge_d, halo_edge_d, disk_edge_d, bulge_edge_d, halo_const_d;
__constant__ double Q_total_d, X_total_d;
__constant__ double sin_inclination_d, cos_inclination_d;

//Galaxy parameters
__constant__ double c_halo_d;
__constant__ double v_halo_d;
__constant__ double a_halo_d;
__constant__ double drtrunc_halo_d;
__constant__ double cusp_d;
__constant__ double halo_stream_d;
__constant__ double n_sersic_d;
__constant__ double ppp_d;
__constant__ double b_n_d;
__constant__ double v_bulge_d;
__constant__ double a_bulge_d;
__constant__ double bulge_stream_d;
__constant__ double rho_0_d;
__constant__ double bh_mass_d;
__constant__ double Plcon_d[100];
__constant__ double M_Disk[DISK_LIMIT];
__constant__ double R_Disk[DISK_LIMIT];
__constant__ double Z_Disk[DISK_LIMIT];
__constant__ double Out_Disk[DISK_LIMIT];
__constant__ double Dr_Trunc[DISK_LIMIT];
__constant__ double R_Kormendy[DISK_LIMIT];
__constant__ double Alpha[DISK_LIMIT];
__constant__ double Sigma_0[DISK_LIMIT];
__constant__ double R_Sigma[DISK_LIMIT];
__constant__ double M_GasDisk[DISK_LIMIT];
__constant__ double R_GasDisk[DISK_LIMIT];
__constant__ double Out_GasDisk[DISK_LIMIT];
__constant__ double Z_GasDisk[DISK_LIMIT];
__constant__ double Dr_Trunc_Gas[DISK_LIMIT];
__constant__ double Gamma[DISK_LIMIT];
__constant__ double R_Kormendy_Gas[DISK_LIMIT];
__constant__ double Alpha_Gas[DISK_LIMIT];
__constant__ double Sigma_0_Gas[DISK_LIMIT];
__constant__ double R_Sigma_Gas[DISK_LIMIT];
__constant__ double Rho_Disk_Const_d[DISK_LIMIT];
__constant__ double Disk_Const_d[DISK_LIMIT];
__constant__ double Rho_GasDisk_Const_d[DISK_LIMIT];
__constant__ double GasDisk_Const_d[DISK_LIMIT];

//Astronomical parameters
__constant__ double ML_Disk[DISK_LIMIT];
__constant__ double ml_bulge_d;
__constant__ double inclination_d;
__constant__ double distance_d;
__constant__ double v_sys_d;
__constant__ double r_sys_d;

// __constant__ thrust::device_ptr<double> Radius_CUDA;
// __constant__ thrust::device_ptr<double> DF_Sersic_CUDA;
// __constant__ thrust::device_ptr<double> DF_NFW_CUDA;
// __constant__ thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
// __constant__ thrust::device_ptr<double> Dens_Psi_Halo_CUDA;
// __constant__ thrust::device_ptr<double> Omega_CUDA;
// __constant__ thrust::device_ptr<double> Omega2_CUDA;
// __constant__ thrust::device_ptr<double> A_K_CUDA;
// __constant__ thrust::device_ptr<double> A_K2_CUDA;

__constant__ int l_d;

//Device function declarations
__host__ __device__ __noinline__ 
double PotCUDA(double s, double z,
               thrust::device_ptr<double> Radius_CUDA,
               thrust::device_ptr<double> A_Pot_CUDA);
__host__ __device__ __noinline__ 
double LegendreCUDA(const int l, const double x);
__host__ __device__ __noinline__ 
double GetTruncCUDA(const double &rad, const int &i);
__host__ __device__ __noinline__     
double GetTruncGasCUDA(const double &rad, const int &i);
__host__ __device__ __noinline__ 
void   GetTruncPrimeCUDA(const double &r, const int &i, 
                         double &truncfac, double &truncfacprime);
__host__ __device__ __noinline__     
void GetTruncGasPrimeCUDA(const double &r, const int &i, 
                          double &truncfac, double &truncfacprime);
__host__ __device__ __noinline__ 
void   DiskProfileCUDA(const double radius, const double z, const int i, 
                       double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__ 
void   DiskProfilePrimeCUDA(const double &radius, const double &z, const int &i, 
                            double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__ 
void   DiskProfile2PrimeCUDA(const double &radius, const double &z, const int &i, 
                             double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__ 
double ForceCUDA(double s, double z, double &force_s, double &force_z,
                 thrust::device_ptr<double> A_Pot_CUDA, 
                 thrust::device_ptr<double> F_R_CUDA);
__host__ __device__ __noinline__ 
void   AppDiskForceCUDA(double &r, double &z, double &fsad, double &fzad);
__host__ __device__ __noinline__ 
double AppDiskPotCUDA(double &r, double &z);
__host__ __device__ __noinline__ 
double AppDiskDensCUDA(double &r, double &z);

__host__ __device__ __noinline__ 
void   ForceCUDA(double s, double z, double &force_s, double &force_z,
                 thrust::device_ptr<double> Radius_CUDA,
                 thrust::device_ptr<double> A_Pot_CUDA, 
                 thrust::device_ptr<double> F_R_CUDA);
__host__ __device__ __noinline__ 
void RotateCoordinatesCUDA(double &xp, double &yp, double &zp, 
                           double &x, double &y, double &z);
__host__ __device__ __noinline__ 
void RotateCoordinatesBackCUDA(double &xp, double &yp, double &zp, 
                               double &x, double &y, double &z);
__host__ __device__ __noinline__ 
void RotateCoordinatesOffCentreCUDA(double &xp, double &yp, double &zp, 
                                    double &x, double &y, double &z);
__host__ __device__ __noinline__ 
double DiskSurfaceDensityCUDA(double &xp, double &yp, int &j,
                              thrust::device_ptr<double> Radius_CUDA, 
                              thrust::device_ptr<double> A_Pot_CUDA);
__host__ __device__ __noinline__ 
double BulgeSurfaceDensityCUDA(double &xp, double &yp, 
                               thrust::device_ptr<double> Radius_CUDA, 
                               thrust::device_ptr<double> A_Pot_CUDA,
                               thrust::device_ptr<double> Dens_Psi_Bulge_CUDA);
__host__ __device__ __noinline__ 
void PreDiskVelocitiesCUDA(double *pos, double &v_phi_max, double *v_cyl, 
                           double &f_max, double &v_circ, int &j, 
                           thrust::device_ptr<double> Radius_CUDA,
                           thrust::device_ptr<double> Rad_Spline_CUDA,
                           thrust::device_ptr<double> A_Pot_CUDA,
                           thrust::device_ptr<double> F_R_CUDA,
                           thrust::device_ptr<double> Omega_CUDA,
                           thrust::device_ptr<double> Omega2_CUDA,
                           thrust::device_ptr<double> A_K_CUDA,
                           thrust::device_ptr<double> A_K2_CUDA,
                           thrust::device_ptr<double> Am_Tab_CUDA,
                           thrust::device_ptr<double> R_Tab_CUDA,
                           thrust::device_ptr<double> R_Tab2_CUDA,
                           thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                           thrust::device_ptr<double> FD_Rat_CUDA,
                           thrust::device_ptr<double> D_Rat2_CUDA,
                           thrust::device_ptr<double> FSZ_Rat_CUDA,
                           thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
void DiskVelocitiesCUDA(double *pos, double &v_phi_max, double *v_cyl, double &f_max,
                        double &v_circ, double &w_los, double &df, int &j, 
                        thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> Rad_Spline_CUDA,
                        thrust::device_ptr<double> A_Pot_CUDA,
                        thrust::device_ptr<double> Omega_CUDA,
                        thrust::device_ptr<double> Omega2_CUDA,
                        thrust::device_ptr<double> A_K_CUDA,
                        thrust::device_ptr<double> A_K2_CUDA,
                        thrust::device_ptr<double> Am_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab2_CUDA,
                        thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                        thrust::device_ptr<double> FD_Rat_CUDA,
                        thrust::device_ptr<double> D_Rat2_CUDA,
                        thrust::device_ptr<double> FSZ_Rat_CUDA,
                        thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
void DiskDispersionCUDA(double &xp, double &yp, double &wbar_disk, 
                        double &disk_dispersion, int &j,
                        thrust::device_ptr<double> Radius_CUDA,
                        thrust::device_ptr<double> Rad_Spline_CUDA,
                        thrust::device_ptr<double> A_Pot_CUDA,
                        thrust::device_ptr<double> F_R_CUDA,
                        thrust::device_ptr<double> Omega_CUDA,
                        thrust::device_ptr<double> Omega2_CUDA,
                        thrust::device_ptr<double> A_K_CUDA,
                        thrust::device_ptr<double> A_K2_CUDA,
                        thrust::device_ptr<double> Am_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab2_CUDA,
                        thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                        thrust::device_ptr<double> FD_Rat_CUDA,
                        thrust::device_ptr<double> D_Rat2_CUDA,
                        thrust::device_ptr<double> FSZ_Rat_CUDA,
                        thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
void BulgeVelocitiesCUDA(double &psi, double &vmag, double &dvmag, double &x, 
                         double &y, double &z, double &w_los, double &dfn,
                         thrust::device_ptr<double> DF_Sersic_CUDA);
__host__ __device__ __noinline__ 
void BulgeDispersionCUDA(double &xp, double &yp, double &wbar_bulge, 
                         double &bulge_dispersion,
                         const thrust::device_ptr<double> Radius_CUDA,
                         const thrust::device_ptr<double> A_Pot_CUDA, 
                         const thrust::device_ptr<double> Dens_Psi_Bulge_CUDA,
                         const thrust::device_ptr<double> DF_Sersic_CUDA);
__host__ __device__ __noinline__ 
void FindMaxCUDA(double &r, double &z, double &vsigp, double &vpmax, 
                 double &fmax, int &j,
                 thrust::device_ptr<double> Radius_CUDA,
                 thrust::device_ptr<double> Rad_Spline_CUDA,
                 thrust::device_ptr<double> A_Pot_CUDA,
                 thrust::device_ptr<double> Omega_CUDA,
                 thrust::device_ptr<double> Omega2_CUDA,
                 thrust::device_ptr<double> A_K_CUDA,
                 thrust::device_ptr<double> A_K2_CUDA,
                 thrust::device_ptr<double> Am_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab2_CUDA,
                 thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                 thrust::device_ptr<double> FD_Rat_CUDA,
                 thrust::device_ptr<double> D_Rat2_CUDA,
                 thrust::device_ptr<double> FSZ_Rat_CUDA,
                 thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
void SplintDCUDA(thrust::device_ptr<double> xa, thrust::device_ptr<double> ya, 
                 thrust::device_ptr<double> y2a, int n, double x, double &y);
__host__ __device__ __noinline__ 
double DiskDensfICUDA(double r, double z, int i, thrust::device_ptr<double> Radius_CUDA, 
                                                 thrust::device_ptr<double> A_Pot_CUDA);
__host__ __device__ __noinline__ 
double BulgeDFCUDA(double &energy, thrust::device_ptr<double> DF_Sersic_CUDA);
__host__ __device__ __noinline__ 
void GetOmegaKappaCUDA(double &r, double &freq_omega, double &freq_kappa,
                       thrust::device_ptr<double> Radius_CUDA,
                       thrust::device_ptr<double> Omega_CUDA,
                       thrust::device_ptr<double> Omega2_CUDA,
                       thrust::device_ptr<double> A_K_CUDA,
                       thrust::device_ptr<double> A_K2_CUDA);
__host__ __device__ __noinline__ 
double BulgeDensPsiCUDA(double energy, thrust::device_ptr<double> Dens_Psi_Bulge_CUDA);
__host__ __device__ __noinline__ 
double DiskDF5ezCUDA(double &vr, double &vt, double &vz, double &r, double &z, int &i,
                     thrust::device_ptr<double> Radius_CUDA,
                     thrust::device_ptr<double> Rad_Spline_CUDA,
                     thrust::device_ptr<double> A_Pot_CUDA,
                     thrust::device_ptr<double> Omega_CUDA,
                     thrust::device_ptr<double> Omega2_CUDA,
                     thrust::device_ptr<double> A_K_CUDA,
                     thrust::device_ptr<double> A_K2_CUDA,
                     thrust::device_ptr<double> Am_Tab_CUDA,
                     thrust::device_ptr<double> R_Tab_CUDA,
                     thrust::device_ptr<double> R_Tab2_CUDA,
                     thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                     thrust::device_ptr<double> FD_Rat_CUDA,
                     thrust::device_ptr<double> D_Rat2_CUDA,
                     thrust::device_ptr<double> FSZ_Rat_CUDA,
                     thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
double DiskDF3ezCUDA(double &ep, double &am, double &ez, int &i,
                     thrust::device_ptr<double> Radius_CUDA,
                     thrust::device_ptr<double> Rad_Spline_CUDA,
                     thrust::device_ptr<double> A_Pot_CUDA,
                     thrust::device_ptr<double> Omega_CUDA,
                     thrust::device_ptr<double> Omega2_CUDA,
                     thrust::device_ptr<double> A_K_CUDA,
                     thrust::device_ptr<double> A_K2_CUDA,
                     thrust::device_ptr<double> Am_Tab_CUDA,
                     thrust::device_ptr<double> R_Tab_CUDA,
                     thrust::device_ptr<double> R_Tab2_CUDA,
                     thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                     thrust::device_ptr<double> FD_Rat_CUDA,
                     thrust::device_ptr<double> D_Rat2_CUDA,
                     thrust::device_ptr<double> FSZ_Rat_CUDA,
                     thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
double FnaMidDenCUDA(double &r, int &j, thrust::device_ptr<double> Radius_CUDA, 
                     thrust::device_ptr<double> Rad_Spline_CUDA,
                     thrust::device_ptr<double> A_Pot_CUDA, 
                     thrust::device_ptr<double> FD_Rat_CUDA,
                     thrust::device_ptr<double> D_Rat2_CUDA);
__host__ __device__ __noinline__ 
double RCircCUDA(double &am, thrust::device_ptr<double> Am_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab_CUDA, 
                 thrust::device_ptr<double> R_Tab2_Zero_CUDA);
__host__ __device__ __noinline__ 
double SigZ2CUDA(double r, int &i, thrust::device_ptr<double> Radius_CUDA, 
                 thrust::device_ptr<double> Rad_Spline_CUDA,
                 thrust::device_ptr<double> A_Pot_CUDA, 
                 thrust::device_ptr<double> FSZ_Rat_CUDA,
                 thrust::device_ptr<double> SZ_Rat2_CUDA);
__host__ __device__ __noinline__ 
double SigR2CUDA(double r, int &i);

__host__ __device__ __noinline__
double GetNFWd2rhodpsi2CUDA(double &r, 
                            thrust::device_ptr<double> Radius_CUDA,
                            thrust::device_ptr<double> H_FR_CUDA,
                            thrust::device_ptr<double> B_FR_CUDA,
                            thrust::device_ptr<double> D_FR_CUDA,
                            thrust::device_ptr<double> D_Dens_CUDA,
                            thrust::device_ptr<double> Halo_AC_Radius_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_D_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA);
__host__ __device__ __noinline__
double GetSersicd2rhodpsi2CUDA(double &r, 
                               thrust::device_ptr<double> Radius_CUDA,
                               thrust::device_ptr<double> H_FR_CUDA,
                               thrust::device_ptr<double> B_FR_CUDA,
                               thrust::device_ptr<double> D_FR_CUDA,
                               thrust::device_ptr<double> D_Dens_CUDA,
                               thrust::device_ptr<double> Halo_AC_Radius_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_D_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA);
__host__ __device__ __noinline__
void FindBracketsCUDA(double &psi, double &r_lower, double &r_upper, 
                      thrust::device_ptr<double> Radius_CUDA,
                      thrust::device_ptr<double> H_Pot_CUDA,
                      thrust::device_ptr<double> B_Pot_CUDA,
                      thrust::device_ptr<double> D_Pot_CUDA);
__host__ __device__ __noinline__
double RootBisectionCUDA(double &psi, double &r_upper, double &r_lower, 
                         double tolerance_d, thrust::device_ptr<double> Radius_CUDA,
                         thrust::device_ptr<double> H_Pot_CUDA,
                         thrust::device_ptr<double> B_Pot_CUDA,
                         thrust::device_ptr<double> D_Pot_CUDA);
__host__ __device__ __noinline__
double GetTotalPsiCUDA(double &r, thrust::device_ptr<double> Radius_CUDA,
                       thrust::device_ptr<double> H_Pot_CUDA,
                       thrust::device_ptr<double> B_Pot_CUDA,
                       thrust::device_ptr<double> D_Pot_CUDA);
__host__ __device__ __noinline__
double DiskForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> D_FR_CUDA);
__host__ __device__ __noinline__
double HaloDensityCUDA(double &r, 
                       thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                       thrust::device_ptr<double> Halo_AC_Dens_CUDA);
__host__ __device__ __noinline__
double HaloDensityPrimeCUDA(double &r, 
                            thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                            thrust::device_ptr<double> Halo_AC_Dens_CUDA, 
                            thrust::device_ptr<double> Halo_AC_Dens_D_CUDA);
__host__ __device__ __noinline__
double HaloDensity2PrimeCUDA(double &r, 
                             thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_D_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA);
__host__ __device__ __noinline__
double RawHaloProfileCUDA(double &r);
__host__ __device__ __noinline__
double RawHaloProfilePrimeCUDA(double &r);
__host__ __device__ __noinline__
double RawHaloProfile2PrimeCUDA(double &r);
__host__ __device__ __noinline__
double HaloProfileDensCUDA(double &r, 
                           thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                           thrust::device_ptr<double> Halo_AC_Dens_CUDA);
__host__ __device__ __noinline__
double HaloProfileDensPrimeCUDA(double &r, 
                                thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                                thrust::device_ptr<double> Halo_AC_Dens_D_CUDA);
__host__ __device__ __noinline__
double HaloProfileDens2PrimeCUDA(double &r, 
                                 thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                                 thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA);
__host__ __device__ __noinline__
double GetHaloTruncCUDA(double &r);
__host__ __device__ __noinline__
double GetHaloTruncPrimeCUDA(double &r);
__host__ __device__ __noinline__
double GetHaloTrunc2PrimeCUDA(double &r);
__host__ __device__ __noinline__
double HaloForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> H_FR_CUDA);
__host__ __device__ __noinline__
double SersicPotentialCUDA(double &r);
__host__ __device__ __noinline__
double SersicMassCUDA(double &r);
__host__ __device__ __noinline__
double SersicDensProfileCUDA(double &r, double *profile);
__host__ __device__ __noinline__
double SersicDensCUDA(double &r);
__host__ __device__ __noinline__
double SersicDensPrimeCUDA(double &r);
__host__ __device__ __noinline__
double SersicDens2PrimeCUDA(double &r);
__host__ __device__ __noinline__
double SersicForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA,
                                  thrust::device_ptr<double>B_FR_CUDA);
__host__ __device__ __noinline__
double DiskForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                     thrust::device_ptr<double> D_Dens_CUDA);
__host__ __device__ __noinline__
double DiskDensityCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                       thrust::device_ptr<double> D_Dens_CUDA);
__host__ __device__ __noinline__     
double GasDiskSurfaceDensfICUDA(double r, int j, 
                                thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> A_Pot_CUDA, 
                                thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double GasDiskSurfaceDensfICUDA2(double r, int j, double dens_const,
                                 thrust::device_ptr<double> Radius_CUDA, 
                                 thrust::device_ptr<double> A_Pot_CUDA, 
                                 thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                 thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double IterateGasDensCUDA(double r, int j,
                          thrust::device_ptr<double> Radius_CUDA, 
                          thrust::device_ptr<double> A_Pot_CUDA, 
                          thrust::device_ptr<double> Polytrope_Const_CUDA, 
                          thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
void GasDiskProfileCUDA(const double radius, const double z, const int i, 
                        double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__     
void GasDiskProfilePrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__     
void GasDiskProfile2PrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile);
__host__ __device__ __noinline__     
void AppGasDiskForceCUDA(double &r, double &z, double &fsad, double &fzad);
__host__ __device__ __noinline__     
double AppGasDiskPotCUDA(double &r, double &z);
__host__ __device__ __noinline__     
double AppGasDiskDensCUDA(double &r, double &z);
__host__ __device__ __noinline__     
double GasDiskDensfICUDA(double r, double z, int i, 
                         thrust::device_ptr<double> Radius_CUDA, 
                         thrust::device_ptr<double> A_Pot_CUDA, 
                         thrust::device_ptr<double> Polytrope_Const_CUDA, 
                         thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double GasDiskDensfINoTruncCUDA(double r, double z, int i, 
                                thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> A_Pot_CUDA, 
                                thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double GasDiskDensfINoTruncCUDA2(double r, double z, int i, double dens_const,
                                 thrust::device_ptr<double> Radius_CUDA, 
                                 thrust::device_ptr<double> A_Pot_CUDA, 
                                 thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                 thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double GetDensConstCUDA(double r, int j, thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> GasDensity_Const_CUDA);
__host__ __device__ __noinline__     
double GetPolyConstCUDA(double r, int j, thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> GasDensity_Const_CUDA);
#endif
