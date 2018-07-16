#include "galaxy.h"

double HaloDens(double r, double z)
{
    double psi = Pot(r,z);
    return HaloDensPsi(psi);
}

double HaloDensPsi(double &energy)
{
    if (energy < psi_crit) return 0;
    
    if (energy >= psi_0) 
    {
        return Dens_Psi_Halo[0];
    }
    
    double rj1 = (psi_0 - energy)/(psi_0 - psi_d);
    //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
    double rj = (n_psi-1.0)*log(rj1)*log_rj2;// /log(rj2);
    //if(nbody_flag==1)
    //cout<<"   "<<rj1<<" "<<rj2<<" "<<rj<<" "<<energy<<" "<<psi_0<<" "<<psi_d<<" "<<psi_crit<<endl;
    int j = int(rj);
    
    if (j < 0)
        j = 0;
    else if (j >= n_psi-1)
        j = n_psi - 2;
    
    double frac = rj - j;
    
    double halo_dens_psi = Dens_Psi_Halo[j] + frac*(Dens_Psi_Halo[j+1]-Dens_Psi_Halo[j]);
    //if(nbody_flag==1)
    //cout<<"halodenspsi  "<<energy<<" "<<j<<" "<<frac<<" "<<halo_dens_psi<<" "<<Dens_Psi_Halo[j+1]<<" "<<Dens_Psi_Halo[j]<<endl;
    
    return halo_dens_psi;
}

