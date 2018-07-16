#ifndef _chisq_h_
#define _chisq_h_

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
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#define data_sets 4
#define m_scale 2.325e9
#define v_scale 100.0
#define mag_sun 4.08

using namespace::std;

extern int do_chisq_file_io, n_los, n_vel;

extern vector<double> Error_Factors;
extern vector<double> Chi_Square;
extern vector<int>    DOF;

extern vector<double> SB_Radii_Orig, P_Angle_Orig, Ellipticity_Orig, SB_Data_Orig, SB_Error_Orig;
extern vector<double> VC_Radii_Orig, VC_Data_Orig, VC_Error_Orig;
extern vector<double> SVel_Radii_Orig, SVel_Data_Orig, SVel_Error_Orig;
extern vector<double> SDisp_Radii_Orig, SDisp_Data_Orig, SDisp_Error_Orig;

extern vector<double> SB_Radii, P_Angle, Ellipticity, SB_Data, SB_Error;
extern vector<double> VC_Radii, VC_Data, VC_Error;
extern vector<double> SVel_Radii, SVel_Data, SVel_Error;
extern vector<double> SDisp_Radii, SDisp_Data, SDisp_Error;

double ArcsecToKpc(double &rad);
void   BulgeDispersion(double &xp, double &yp, double &wbar_bulge, 
                       double &bulge_dispersion);
double BulgeSurfaceDensity(double &xp, double &yp);
void   BulgeVelocities(double &psi, double &vmag, double &dvmag, double &x, 
                       double &y, double &z, double &w_los, double &dfn);
void   ConvertRadius(vector<double> &Vec_Orig, vector<double> &Vec);
void   ConvertVelocities(vector<double> &Vec1_Orig, vector<double> &Vec2_Orig, 
                         vector<double> &Vec1, vector<double> &Vec2);
void   DiskDispersion(double &xp, double &yp, double &wbar_disk, 
                      double &disk_dispersion, int &j);
double DiskSurfaceDensity(double &xp, double &yp, int &j);
void   DiskVelocities(double *pos, double &v_phi_max, double *v_cyl, double &f_max, 
                      double &v_circ, double &w_los, double &df, int &j);
void   FindMax(double &r, double &z, double &vsigp, double &vpmax, 
               double &fmax, int &j);
void   FindMaxGas(double &r, double &z, double &vsigp, double &vpmax, 
                  double &fmax, int &j);
void   GasDiskDispersion(double &xp, double &yp, double &wbar_disk, 
                         double &disk_dispersion, int &j);
double GasDiskSurfaceDensity(double &xp, double &yp, int &j);
void   GasDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, double &f_max,
                         double &v_circ, double &w_los, double &df, int &j);
double GetChiSquare(vector<double> &Chi_Square);
double GetCircVelocity(vector<double>& Radii, vector<double>& Data, 
                       vector<double>& Error);
double GetHIDispersion(vector<double>& Radii, vector<double>& Data, 
                       vector<double>& Error);
double GetHIVelocity(vector<double>& Radii, vector<double>& Data,
                     vector<double>& Error);
double GetStarDispersion(vector<double>& Radii, vector<double>& Data, 
                         vector<double>& Error);
double GetStarVelocity(vector<double>& Radii, vector<double>& Data, 
                       vector<double>& Error);
double GetSurfaceBrightness(vector<double>& Radii, vector<double>& Data, 
                            vector<double>& Error);
void   GetObservedData(void);
void   PreDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, double &f_max, 
                         double &v_circ, int &j);
void   PreGasDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, 
                            double &f_max, double &v_circ, int &j);
void   ReadData(ifstream& file, vector<double>& Vec_1, 
                vector<double>& Vec_2, vector<double>& Vec_3);
void   ReadData(ifstream& file, vector<double>& Vec_1, vector<double>& Vec_2, 
                vector<double>& Vec_3, vector<double>& Vec_4);
void   ReadData(ifstream& file, vector<double>& Vec_1, vector<double>& Vec_2, 
                vector<double>& Vec_3, vector<double>& Vec_4, vector<double>& Vec_5);
void   RotateCoordinates(double &xp, double &yp, double &zp, double &x, 
                         double &y, double &z);
void   RotateCoordinatesBack(double &xp, double &yp, double &zp, double &x, 
                             double &y, double &z);

#endif
