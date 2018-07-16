#include "galaxy.h"

double SigZ2(double r, int &i)
{
    double psi_zh = Pot(r, G.Z_Disk[i]);
    double psi_0 = Pot(r, 0.0);
    double true_sig_z2 = (psi_zh-psi_0)/log(0.419974341614);
    double f_cor;
    
    SplintD(Rad_Spline[i], FSZ_Rat[i], SZ_Rat2[i], nr_spline, r, f_cor);
    //cout << "sigz   " << FSZ_Rat[0] << " " << FSZ_Rat[1] << " " << SZ_Rat2[0] 
    //     << " " << SZ_Rat2[1] << endl;
    
    return f_cor*true_sig_z2;
}

double SigR2(double r, int &i)
{
    return G.Sigma_0[i]*G.Sigma_0[i]*exp(-r/G.R_Sigma[i]);
}

double SigZ2Gas(double r, int &i)
{
//     double psi_zh = Pot(r, G.Z_GasDisk[i]);
//     double psi_0 = Pot(r, 0.0);
//     double true_sig_z2 = (psi_zh-psi_0)/log(0.419974341614);
//     double f_cor;
//     
//     SplintD(Rad_Spline_Gas[i], FSZ_Rat_Gas[i], SZ_Rat2_Gas[i], nr_spline, r, f_cor);
//     //cout << "sigz   " << FSZ_Rat[0] << " " << FSZ_Rat[1] << " " 
       //     << SZ_Rat2[0] << " " << SZ_Rat2[1] << endl;
//     return f_cor*true_sig_z2;
    
    return G.Sigma_0_Gas[i]*G.Sigma_0_Gas[i];
    
    //SplintD(Rad_Spline_Gas[i], FSZ_Rat_Gas[i], SZ_Rat2_Gas[i], nr_spline, r, f_cor);
    //cout << "sigz   " << FSZ_Rat[0] << " " << FSZ_Rat[1] << " " << SZ_Rat2[0] << " " << SZ_Rat2[1] << endl;
    //return f_cor*true_sig_z2;
}

double SigR2Gas(double r, int &i)
{
    return G.Sigma_0_Gas[i]*G.Sigma_0_Gas[i]*exp(-r/G.R_Sigma_Gas[i]);
}

//The totals are the mass weighted sums of the component dispersions because
//<sigma>=0 but <sigma^2>!=0
double SigZ2Total(double r)
{
    double total_sig = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double sig_z = SigZ2(r,i);
        total_sig += G.M_Disk[i]/total_disk_mass*sig_z*sig_z;
    }
    
    return sqrt(total_sig);
}

double SigR2Total(double r)
{
    double total_sig = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double sig_r = SigR2(r,i);
        total_sig += G.M_Disk[i]/total_disk_mass*sig_r*sig_r;
    }
    
    return sqrt(total_sig);
}
