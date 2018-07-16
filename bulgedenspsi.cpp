#include "galaxy.h"

double BulgeDens(double r, double z)
{
    double psi = Pot(r, z);
    return BulgeDensPsi(psi);
}

double BulgeDensPsi(double &energy)
{
    if (energy < psi_crit) return 0;
    
    else if (energy >= psi_0) return Dens_Psi_Bulge[0];
    
    double rj1 = (psi_0 - energy)/(psi_0 - psi_d);
    //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
    double rj = (n_psi-1.0)*log(rj1)*log_rj2;
    
	int j = int(rj);
    
    if (j < 0)
        j = 0;
    else if (j >= n_psi-1)
        j = n_psi - 2;
    
    double frac = rj - j;
    
    //cout << "bdpsi " << rj1 << " " << rj2 << " " << rj << " " << j << 
    //        " " << setprecision(12) << psi_0 << " " << energy << endl;
    
    return Dens_Psi_Bulge[j] + frac*(Dens_Psi_Bulge[j+1] - Dens_Psi_Bulge[j]);
}
    
    
