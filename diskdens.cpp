//Get the disk density from the radius and height above the plane

#include "galaxy.h"

double DiskDens(double r, double z, double psi)
{
    double diskdens = 0;
    
    if (!disk_flag)
    {
        return diskdens = 0;
    }
    
    for (int i = 0; i < disk; ++i)
    {
        if(fabs(z/G.Z_Disk[i])>30)
            continue;
        
        double trunc_fac = GetTrunc(r, i);
        
        if (trunc_fac==0)
            continue;

        double con;
        
        if (z==0)
        {
            con = 1;
        }
        else if (!halo_flag && !bulge_flag)
        {
            con = exp(-fabs(z/G.Z_Disk[i]));
            con = 2.0*con/(1.0+con*con);
            con *= con;
        }
        else
        {
            //KK: make disk vertically isothermal; normalize density so that
            //at height 3*zdisk the vertical density drop is that of a sech^2.
            double psizh = Pot(r, 3.0*G.Z_Disk[i]);
            double psi00 = Pot(r, 0);
            double dpsizh = psi00 - psizh;
            double dpsi = psi00 - psi;
            double coeff = 0;
            
            if (r==0 && z==0)
            {
                dpsi = 0;
            }

            if (fabs(dpsizh) > 0)
            {
                coeff = dpsi / dpsizh;
            }
            
            if (coeff > 16 || coeff < 0)
            {
                con = 0;
            }
            else
            {
                //the 0.009866 is from changing from Z_Disk to 3*Z_Disk
                //Otherwise it would exp(-0.8676)=0.419974
                con = pow(0.009866, coeff);
            }
        }

		double Sigma_Profile[2], Rho_Profile[2];
		
        DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        diskdens += Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
//         double zmaximum = 2*G.Z_Disk[i];
//         double facfac = zmaximum*(tanh((r-3*G.R_Disk[i])/G.R_Disk[i])+1.0);
//         diskdens += Rho_Disk_Const[i]*Sigma_Profile[0]*trunc_fac/
//                     (cosh(z/facfac)*cosh(z/facfac));
    }
    //if(r>299||z>299)cout<<"diskdens   "<<r<<" "<<z<<" "<<diskdens<<endl;
    return diskdens;
}

double DiskDensf(double r, double z)
{
    double psi = Pot(r, z);
    return DiskDens(r, z, psi);
}

double DiskDensPsi(double r, double z, double psi)
{
    double diskdenspsi = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double con;
        double trunc_fac = GetTrunc(r, i);
    
        if (z == 0)
            con = 1;
        else
        {
            double psizh = Pot(r, G.Z_Disk[i]);
            double psi00 = Pot(r, 0);
            con = pow(0.419974341614, (psi-psi00)/(psizh-psi00));
            //con = exp((psi-psi00)/0.0144);
            //con *= (1-exp(-r*r/5));
        
        }
            
		double Sigma_Profile[2], Rho_Profile[2];
		
        DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
        
        diskdenspsi += Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
//         double zmaximum = 2*G.Z_Disk[i];
//         double facfac = zmaximum*(tanh((r-3*G.R_Disk[i])/G.R_Disk[i])+1.0);
//         diskdenspsi += Rho_Disk_Const[i]*Sigma_Profile[0]*trunc_fac/
//                     (cosh(z/facfac)*cosh(z/facfac));
    }
    
    return diskdenspsi;
}

double DiskDensI(double r, double z, double psi, int i)
{
    double diskdens;
    
    if (!disk_flag)
    {
        return diskdens = 0;
    }
    
    double con;
    double trunc_fac = GetTrunc(r, i);

    if (!halo_flag && !bulge_flag && !gasdisk_flag)
    {
        con = exp(-fabs(z/G.Z_Disk[i]));
        con = pow(2.0*con / (1.0+con*con), 2.0);
    }
    else
    {
        double psizh = Pot(r, 3.0*G.Z_Disk[i]);
        double psi00 = Pot(r, 0);
        double dpsizh = psi00 - psizh;
        double dpsi = psi00 - psi;
        double coeff = 0;

        if (r==0 && z==0)
        {
            dpsi = 0;
        }

        if (fabs(dpsizh) > 0)
        {
            coeff = dpsi / dpsizh;
        }

        if (coeff > 16 || coeff < 0)
        {
            con = 0;
        }
        else
        {
            con = pow(0.009866, coeff);
        }
            
            //con *= (1-exp(-r*r/5));
            //con = exp(-dpsi/0.0144);
    }

	double Sigma_Profile[2], Rho_Profile[2];

    DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

    //double zfac = G.Z_Disk[i]*(1-0.95*exp(-r*r/25));//+0.00001;
    //con = 1.0/cosh(fabs(z/zfac))/cosh(fabs(z/zfac));
    diskdens = Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
//         double zmaximum = 2*G.Z_Disk[i];
//         double facfac = zmaximum*(tanh((r-3*G.R_Disk[i])/G.R_Disk[i])+1.0);
//         diskdens += Rho_Disk_Const[i]*Sigma_Profile[0]*trunc_fac/
//                     (cosh(z/facfac)*cosh(z/facfac));
    
    return diskdens;
}

double DiskDensfI(double r, double z, int i)
{
    double psi = Pot(r, z);
    return DiskDensI(r, z, psi, i);
}

//Now do the gas density

double GasDiskDens(double r, double z, double psi)
{
    double psi_zh = Pot(r, G.Z_Disk[0]);
    double psi_0 = Pot(r, 0.0);
    double true_stellar_sig_z2 = (psi_zh-psi_0)/log(0.419974);
    double gas_sig_z2 = G.Z_GasDisk[0];
    double diskdens = 0;
    
    double p = true_stellar_sig_z2/gas_sig_z2;
    
    if (!gasdisk_flag)
    {
        return diskdens = 0;
    }
    
    for (int i = 0; i < gas_disk; ++i)
    {
        if(fabs(z/G.Z_GasDisk[i])>30)
            continue;
        
        double trunc_fac = GetTruncGas(r, i);
        
        if (trunc_fac==0)
            continue;

        double con;
        
        if (z==0)
        {
            con = 1;
        }
        else if (!halo_flag && !bulge_flag && !disk_flag)
        {
            con = exp(-fabs(z/G.Z_GasDisk[i]));
            con = 2.0*con/(1.0+con*con);
            con *= con;
        }
        else
        {
            //KK: make disk vertically isothermal; normalize density so that
            //at height 3*zdisk the vertical density drop is that of a sech^2.
            double psizh = Pot(r, 3.0*G.Z_GasDisk[i]);
            double psi00 = Pot(r, 0);
            double dpsizh = psi00 - psizh;
            double dpsi = psi00 - psi;
            double coeff = 0;
            double base = pow(0.009866037, p);//this is (sech^2(3))^p
            
            if (r==0 && z==0)
            {
                dpsi = 0;
            }

            if (fabs(dpsizh) > 0)
            {
                coeff = dpsi / dpsizh;
            }
            
            if (coeff > 16 || coeff < 0)
            {
                con = 0;
            }
            else
            {
                //the 0.009866 is from changing from Z_Disk to 3*Z_Disk
                //Otherwise it would 0.8676
                con = pow(base, coeff);            
            }
        }

		double Sigma_Profile[2], Rho_Profile[2];
		
        DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
        diskdens += Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
    }
    return diskdens;
}

double GasDiskDensf(double r, double z)
{
    double psi = Pot(r, z);
    return GasDiskDens(r, z, psi);
}

double GasDiskDensPsi(double r, double z, double psi)
{
    double psi_zh = Pot(r, G.Z_Disk[0]);
    double psi_0 = Pot(r, 0.0);
    double true_stellar_sig_z2 = (psi_zh-psi_0)/log(0.419974);
    double gas_sig_z2 = G.Z_GasDisk[0];
    
    double p = true_stellar_sig_z2/gas_sig_z2;
    double diskdenspsi = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        double con;
        double trunc_fac = GetTrunc(r, i);
    
        if (z == 0)
            con = 1;
        else
        {
            double psizh = Pot(r, G.Z_GasDisk[i]);
            double psi00 = Pot(r, 0);
            con = pow(0.419974341614, p*(psi-psi00)/(psizh-psi00));
        }
        
		double Sigma_Profile[2], Rho_Profile[2];
		
        GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
        
        diskdenspsi += Rho_GasDisk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
    }
    
    return diskdenspsi;
}

double GasDiskDensI(double r, double z, double psi, int i)
{
    double diskdens;
    
    if (!disk_flag)
    {
        return diskdens = 0;
    }
    
    double con;
    double trunc_fac = GetTrunc(r, i);

    if (!halo_flag && !bulge_flag)
    {
        con = exp(-fabs(z/G.Z_GasDisk[i]));
        con = pow(2.0*con / (1.0+con*con), 2.0);
    }
    else
    {
        double psizh = Pot(r, 3.0*G.Z_GasDisk[i]);
        double psi00 = Pot(r, 0);
        double dpsizh = psi00 - psizh;
        double dpsi = psi00 - psi;
        double coeff = 0;

        if (r==0 && z==0)
        {
            dpsi = 0;
        }

        if (fabs(dpsizh) > 0)
        {
            coeff = dpsi / dpsizh;
        }

        if (coeff > 16 || coeff < 0)
        {
            con = 0;
        }
        else
        {
            con = pow(0.009866, coeff);
        }
    }

	double Sigma_Profile[2], Rho_Profile[2];

    GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

    diskdens = Rho_GasDisk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
    
    return diskdens;
}

double GasDiskDensfI(double r, double z, int i)
{
    double psi = Pot(r, z);
    return GasDiskDensI(r, z, psi, i);
}
