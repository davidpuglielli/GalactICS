//Get the disk and halo density from polar coodinates

#include "galaxy.h"

double dPolarDiskDens(double r, double cos_theta)
{
    double z = r*cos_theta;
    double s = r*sqrt(1.0 - cos_theta*cos_theta);
    
    return DiskDensEstimate(s, z);
}

double DiskDensEstimate(double s, double z)
{
    double r = sqrt(s*s + z*z);
    double dens_estimate = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double trunc_fac = GetTrunc(r, i);
		double Sigma_Profile[2], Rho_Profile[2];
    
        DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
    
        dens_estimate += Rho_Disk_Const[i]*Rho_Profile[0]*trunc_fac;
    }
    
    return dens_estimate;
}

double dPolarGasDiskDens(double r, double cos_theta)
{
    double z = r*cos_theta;
    double s = r*sqrt(1.0 - cos_theta*cos_theta);
    
    return GasDiskDensEstimate(s, z);
}

double GasDiskDensEstimate(double s, double z)
{
    double r = sqrt(s*s + z*z);
    double dens_estimate = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        double trunc_fac = GetTruncGas(r, i);
		double Sigma_Profile[2], Rho_Profile[2];
    
        GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
    
        dens_estimate += Rho_GasDisk_Const[i]*Rho_Profile[0]*trunc_fac;
    }
    
    return dens_estimate;
}

double PolarHaloDens(double r, double cos_theta, int l)
{
    //KK: returns density at r,theta, multiplied by the lth harmonic
    //(must integrate this over theta to get the harmonic coeff)
    double z = r*cos_theta;
    double s = r*sqrt(1 - cos_theta*cos_theta);
    
    return HaloDens(s, z)*Plcon[l]*gsl_sf_legendre_Pl(l, cos_theta);
}

double PolarBulgeDens(double r, double cos_theta, int l)
{
    //KK: returns density at r,theta, multiplied by the lth harmonic
    //(must integrate this over theta to get the harmonic coeff)
    double z = r*cos_theta;
    double s = r*sqrt(1 - cos_theta*cos_theta);
	double pot = Pot(s,z);
    
    return BulgeDensPsi(pot)*Plcon[l]*gsl_sf_legendre_Pl(l, cos_theta);
}

double PolarDens(double r, double cos_theta, int l)
{
    //KK: returns density at r,theta, multiplied by the lth harmonic
    //(must integrate this over theta to get the harmonic coeff)
    double z = r*cos_theta;
    double s = r*sqrt(1-cos_theta*cos_theta);
    double polardens = Dens(s,z);
    double polar_dens = polardens*Plcon[l]*gsl_sf_legendre_Pl(l, cos_theta);
    
    return polar_dens;
}

