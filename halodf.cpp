#include "galaxy.h"

double HaloDF(double energy)
{
    if (energy < psi_crit) return 0;
    
    if (energy >= psi_0) return exp(DF_NFW[0]);
    
    double rj1 = (psi_0 - energy)/(psi_0 - psi_d);
    //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
    double rj = 1 + (n_psi-1.0)*log(rj1)*log_rj2;// /log(rj2);
    
    int j = int(rj);
    
    if (j < 1)
        j = 1;
    else if (j >= n_psi)
        j = n_psi - 1;
    
    double frac = rj - j;
    
    //cout<<"halodf   "<<rj<<" "<<DF_NFW[j-1]<<" "<<DF_NFW[j]<<" "<<frac<<" "<<exp(DF_NFW[j-1] + frac*(DF_NFW[j] - DF_NFW[j-1]))<<endl;
    return exp(DF_NFW[j-1] + frac*(DF_NFW[j] - DF_NFW[j-1]));
}

