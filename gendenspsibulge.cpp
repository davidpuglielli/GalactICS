//Get the bulge density from the distribution function, using the 
//energy as input

#include "galaxy.h"

void GenDensPsiBulge(void)
{
    ofstream denspsifile;
            
    if (do_file_io)
    {
        denspsifile.open("denspsibulge.dat", ios::out);
    }
    
    //double fcut_bulge = 0;
    
    for (int i = 0; i < n_psi; ++i)
    {
        double psi = Table_E[i];
        Dens_Psi_Bulge[i] = GetBulgeDens(psi);
        if (i == n_psi-1)
            Dens_Psi_Bulge[i] = 0;
        
        if (do_file_io)
        {
            denspsifile << setprecision(12) << Table_E[i] << "   " << Dens_Psi_Bulge[i] << endl;
        }
    }
}

double GetBulgeDens(double &psi)
{
    int m = 100;
    double v02 = G.v_bulge * G.v_bulge;
    double v_min = -10;//actually log vmin

    double v_max = log(sqrt(2*(psi-psi_crit)));
    double dlogv = (v_max - v_min) / (m-1.0);
    double logv, vel, energy, df_lowered, rho = 0;
    
    fcut_bulge = BulgeDF(psi_crit);
    fcut_bulge = 0;
    
    for (int i = 0; i < 4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = BulgeDF(energy) - fcut_bulge;
        if (df_lowered > 0)
            rho += fourpi*dlogv*vel*vel*vel*df_lowered*coef(i+1);
    }
    
    for (int i = 4; i < m-4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = BulgeDF(energy) - fcut_bulge;
        rho += fourpi*dlogv*vel*vel*vel*df_lowered;
    }
    
    for (int j = 4; j > 1; --j)
    {
        int jj = m-j+1;
        logv = v_min + (jj-1)*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = BulgeDF(energy) - fcut_bulge;
        rho += fourpi*dlogv*vel*vel*vel*df_lowered*coef(j);
    }
    
    return rho;
}
