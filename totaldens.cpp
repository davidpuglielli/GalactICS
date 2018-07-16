//Get the total density and the total density minus the approximate disk component

#include "galaxy.h"

double TotalDens(double &r, double &z)
{
    double totdens = 0;
    double psi = Pot(r, z);//if(iter>0)
    
    if (disk_flag == 1)
        totdens += DiskDens(r, z, psi);//
    if (gasdisk_flag)
        totdens += GasDiskDens2(r, z, psi);
    if (halo_flag == 1)
        totdens += HaloDensPsi(psi);
    if (bulge_flag == 1)
        totdens += BulgeDensPsi(psi);
    
    return totdens;
}

double Dens(double r, double z)
{
    double dens = 0;
	
	if (disk_flag && gasdisk_flag)
	{
	    double app_dens = AppDiskDens(r,z);
	    double app_dens_gas = AppGasDiskDens(r,z);
		dens = TotalDens(r,z) - app_dens - app_dens_gas;
    } 
    else if (gasdisk_flag)
    {
	    double app_dens = AppGasDiskDens(r,z);
		dens = TotalDens(r,z) - app_dens;
    }
    else if (disk_flag)
    {
	    double app_dens = AppDiskDens(r,z);
		dens = TotalDens(r,z) - app_dens;
    }
    else
	{
	    dens = TotalDens(r,z);
	}
	
	return dens;
}
