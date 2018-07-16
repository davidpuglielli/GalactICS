//Get the halo density from the distribution function, using the 
//energy as input

#include "galaxy.h"

void GenDensPsiHalo(void)
{
    ofstream denspsifile;
    
    if (do_file_io)
    {
        denspsifile.open("denspsihalo.dat", ios::out);
    }
    //double fcut_halo = 0;
    
    for (int i = 0; i < n_psi; ++i)
    {
        double psi = Table_E[i];
        Dens_Psi_Halo[i] = GetHaloDens(psi);
        if (i == n_psi-1)
            Dens_Psi_Halo[i] = 0;
        
        if (do_file_io)
        {
            denspsifile << setprecision(12) << Table_E[i] << "   " << Dens_Psi_Halo[i] << endl;
        }
        
        //cout<<"gendenspsihalo   "<<log(psi)<<" "<<log(Dens_Psi_Halo[i])<<endl;
    }
}

double GetHaloDens(double &psi)
{
    int m = 100;
    double v_min = -10;//actually log vmin

    double v_max = log(sqrt(2*(psi-psi_crit)));//cout<<"vmax  "<<v_max<<endl;
    double dlogv = (v_max - v_min) / (m-1.0);//cout<<"vmax  "<<dlogv<<endl;
    double logv, vel, energy, df_lowered, rho = 0;
    
    fcut_halo = HaloDF(psi_crit);//cout<<"fcuthalo  "<<fcut_halo<<endl;
    fcut_halo = 0;
    
    for (int i = 0; i < 4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = HaloDF(energy) - fcut_halo;
        if (df_lowered > 0)
            rho += fourpi*dlogv*vel*vel*vel*df_lowered*coef(i+1);
        //cout<<"RHO   "<<vel<<" "<<energy<<" "<<rho<<endl;
    }//cout << "rho      " << rho<<endl;
    
    for (int i = 4; i < m-4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = HaloDF(energy) - fcut_halo;
        rho += fourpi*dlogv*vel*vel*vel*df_lowered;
        //cout<<"RHO   "<<vel<<" "<<energy<<" "<<rho<<endl;
    }//cout << "rho      " << rho<<endl;
    
    for (int j = 4; j > 1; --j)
    {
        int jj = m-j+1;
        logv = v_min + (jj-1)*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        df_lowered = HaloDF(energy) - fcut_halo;
        rho += fourpi*dlogv*vel*vel*vel*df_lowered*coef(j);
        //cout<<"RHO   "<<vel<<" "<<energy<<" "<<rho<<endl;
    }//cout << "rho      " << rho << endl;
    
    return rho;
}
            
