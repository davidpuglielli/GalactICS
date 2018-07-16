//Get the *gas* disk density from the radius and height above the plane
//Note that poly_const is obtained by linear interpolation in the 
//GetPolyConst function.

#include "galaxy.h"

double GasDiskDens(double r, double z, double psi)
{
//     double psi_zh = Pot(r, G.Z_Disk[0]);
//     double psi_0 = Pot(r, 0.0);
//     double true_stellar_sig_z2 = (psi_zh-psi_0)/log(0.419974);
//     double gas_sig_z2 = G.Z_GasDisk[0];
//     double diskdens = 0;
//     
//     double p = true_stellar_sig_z2/gas_sig_z2;
//     
//     if (!gasdisk_flag)
//     {
//         return diskdens = 0;
//     }
//     
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         if(fabs(z/G.Z_GasDisk[i])>30)
//             continue;
//         
//         double trunc_fac = GetTruncGas(r, i);
//         
//         if (trunc_fac==0)
//             continue;
// 
//         double con;
//         
//         if (z==0)
//         {
//             con = 1;
//         }
//         else if (!halo_flag && !bulge_flag && !disk_flag)
//         {
//             con = exp(-fabs(z/G.Z_GasDisk[i]));
//             con = 2.0*con/(1.0+con*con);
//             con *= con;
//         }
//         else
//         {
//             //KK: make disk vertically isothermal; normalize density so that
//             //at height 3*zdisk the vertical density drop is that of a sech^2.
//             double psizh = Pot(r, 3.0*G.Z_GasDisk[i]);
//             double psi00 = Pot(r, 0);
//             double dpsizh = psi00 - psizh;
//             double dpsi = psi00 - psi;
//             double coeff = 0;
//             double base = pow(0.009866037, p);//this is (sech^2(3))^p
//             
//             if (r==0 && z==0)
//             {
//                 dpsi = 0;
//             }
// 
//             if (fabs(dpsizh) > 0)
//             {
//                 coeff = dpsi / dpsizh;
//             }
//             
//             if (coeff > 16 || coeff < 0)
//             {
//                 con = 0;
//             }
//             else
//             {
//                 //the 0.009866 is from changing from Z_Disk to 3*Z_Disk
//                 //Otherwise it would 0.8676
//                 con = pow(base, coeff);            
//             }
//         }
// 
// 		double Sigma_Profile[2], Rho_Profile[2];
// 		
//         DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
//         //if(r>299||z>299)cout<<"dens   "<<r<<" "<<z<<" "<<con<<" "<<Rho_Disk_Const[i]*Sigma_Profile[0]<<" "<<con<<" "<<trunc_fac<<endl;
//         diskdens += Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
//     }
//     //if(r>299||z>299)cout<<"diskdens   "<<r<<" "<<z<<" "<<diskdens<<endl;
//     return diskdens;

    //Interpolate the value of the polytropic constant K
    
    //if (psi < 0)
    //    return 0;
    
//     double psi_edge = Pot(r, Radius[nr-1]);
//     psi -= psi_edge;
//     
//     double density = 0;
//     
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double trunc_fac = GetTruncGas(r, i);
//          
//         if (trunc_fac==0)
//             continue;
//         
//         double polytrope_const = GetPolyConst(r, i);
//                      
//         double exponent = 1/(G.Gamma[i]-1);
//         double constant = (1-G.Gamma[i])/(G.Gamma[i]*polytrope_const);
//         
// 		//double Sigma_Profile[2], Rho_Profile[2];
// 		
//         //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
//         
//         //density += pow(-constant*psi, exponent)*Sigma_Profile[0]*trunc_fac;
//         density += pow(-constant*psi, exponent)*trunc_fac;
//         
//         //cout << "Gasd " << r << " " << pow(-constant*psi, exponent) << " " 
//         //     << "       " << psi << " " << constant << " " << exponent << " " 
//         //     << polytrope_const << endl;
//     }
//         
//     return density;   


    double psi_zero = Pot(r, 0);
    psi -= psi_zero;
    
    double density = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        double trunc_fac = GetTruncGas(r, i);
         
        if (trunc_fac==0)
            continue;
        
        //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double poly_const = GetPolyConst(r, i);
        double dens_const = GetDensConst(r, i);
        
        double rho_0 = dens_const;//*Sigma_Profile[0];

        density += rho_0*exp(psi/poly_const)*trunc_fac;//pow(-constant*psi, exponent)*trunc_fac;//*Sigma_Profile[0];
    }
        
    return density;   
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
    double true_stellar_sig_z2 = (psi_zh-psi_0)/log(0.419974341614);
    double gas_sig_z2 = G.Z_GasDisk[0];
    
    double p = true_stellar_sig_z2/gas_sig_z2;
    double diskdenspsi = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        double con;
        double trunc_fac = GetTruncGas(r, i);
    
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

double GasDiskDensI(double r, double z, double psi, int i, double poly_const)
{
//     double diskdens;
//     
//     if (!disk_flag)
//     {
//         return diskdens = 0;
//     }
//     
//     double con;
//     double trunc_fac = GetTrunc(r, i);
// 
//     if (!halo_flag && !bulge_flag)
//     {
//         con = exp(-fabs(z/G.Z_GasDisk[i]));
//         con = pow(2.0*con / (1.0+con*con), 2.0);
//     }
//     else
//     {
//         double psizh = Pot(r, 3.0*G.Z_GasDisk[i]);
//         double psi00 = Pot(r, 0);
//         double dpsizh = psi00 - psizh;
//         double dpsi = psi00 - psi;
//         double coeff = 0;
// 
//         if (r==0 && z==0)
//         {
//             dpsi = 0;
//         }
// 
//         if (fabs(dpsizh) > 0)
//         {
//             coeff = dpsi / dpsizh;
//         }
// 
//         if (coeff > 16 || coeff < 0)
//         {
//             con = 0;
//         }
//         else
//         {
//             con = pow(0.009866, coeff);
//         }
//     }
// 
// 	double Sigma_Profile[2], Rho_Profile[2];
// 
//     GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
// 
//     diskdens = Rho_GasDisk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
//     
//     return diskdens;

    //if (psi < 0)
    //    return 0;
    
//     double psi_edge = Pot(r, Radius[nr-1]);
//     psi -= psi_edge;
//     
//     double density;
//     
//     double trunc_fac = GetTruncGas(r, i);
// 
//     if (trunc_fac==0)
//     {
//         cout << "Trunc zero " << r << endl;
//         return 0;
//     }
// 
//     double exponent = 1/(G.Gamma[i]-1);
//     double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);
// 
// 	double Sigma_Profile[2], Rho_Profile[2];
// 
//     GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
// 
//     density = pow(-constant*psi, exponent)*trunc_fac;
//     
//     //if (density == 0) 
//     //    cout << "dens zero " << r << " " << z << " " << constant << " " << psi 
//     //         << "                  " << pow(-constant*psi, exponent) << " " 
//     //         << Sigma_Profile[0] << " " << trunc_fac << " " << density << endl;
//         
//     return density;   
    
    double psi_zero = Pot(r, 0);
    psi -= psi_zero;
    
	double Sigma_Profile[2], Rho_Profile[2];

    GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

    double density, rho_0 = Rho_GasDisk_Const[i]*Sigma_Profile[0];
   
    double trunc_fac = GetTruncGas(r, i);

    if (trunc_fac==0)
    {
        //cout << "Trunc zero " << r << endl;
        return 0;
    }

    density = rho_0*exp(psi/poly_const);//pow(-constant*psi, exponent)*trunc_fac;//*Sigma_Profile[0];
    
    if (density < 0)
    {
        //cout << "psi " << r << " " << z << " " << psi << endl; 
        return 0;
    }
        
    return density;   
}

double GasDiskDensfI(double r, double z, int i, double poly_const)
{
    double psi = Pot(r, z);
    return GasDiskDensI(r, z, psi, i, poly_const);
}

double GasDiskDensfI(double r, double z, int i)
{
    double psi = Pot(r, z);
    return GasDiskDensI(r, z, psi, i, 100);
}

// double GasDiskDensINoTrunc(double r, double z, double psi, int i, double poly_const)
// {
// //     //if (psi < 0)
// //     //    return 0;
// //     
// //     double psi_edge = Pot(r, Radius[nr-1]);
// //     psi -= psi_edge;
// //     
// //     double density;
// //     
// //     double trunc_fac = GetTruncGas(r, i);
// // 
// //     if (trunc_fac==0)
// //     {
// //         //cout << "Trunc zero " << r << endl;
// //         return 0;
// //     }
// // 
// //     double exponent = 1/(G.Gamma[i]-1);
// //     double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);
// // 
// // 	double Sigma_Profile[2], Rho_Profile[2];
// // 
// //     GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
// // 
// //     density = pow(-constant*psi, exponent)*trunc_fac;//*Sigma_Profile[0];
// //     
// //     //if (density == 0) 
// //     //    cout << "dens zero " << r << " " << z << " " << constant << " " << psi 
// //     //         << "                  " << pow(-constant*psi, exponent) << " " 
// //     //         << Sigma_Profile[0] << " " << trunc_fac << " " << density << endl;
// //     
// //     if (density < 0)
// //     {
// //         //cout << "psi " << r << " " << z << " " << psi << endl; 
// //         return 0;
// //     }
// //         
// //     return density;   
//     //if (psi < 0)
//     //    return 0;
//     
//     double psi_zero = Pot(r, 0);
//     psi -= psi_zero;
//     
// 	double Sigma_Profile[2], Rho_Profile[2];
// 
//     GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
// 
//     double density, rho_0 = Rho_GasDisk_Const[i]*Sigma_Profile[0];
//    
//     double trunc_fac = GetTruncGas(r, i);
// 
//     if (trunc_fac==0)
//     {
//         //cout << "Trunc zero " << r << endl;
//         return 0;
//     }
// 
//     density = rho_0*exp(psi/poly_const);//*trunc_fac;;//pow(-constant*psi, exponent)*trunc_fac;//*Sigma_Profile[0];
//     
//     if (density <= 0)
//     {
//         //cout << "psi " << r << " " << z << " " << psi << " " << poly_const 
//         //     << " " << exp(psi/poly_const) << endl; 
//         return 0;
//     }
//         
//     return density;   
// }

////////////////////////////////////////////////////////////////////////////////

double GasDiskDens2(double r, double z, double psi)
{
    double psi_zero = Pot(r, 0);
    psi -= psi_zero;
    
    double density = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        if (G.Gamma[i] == 1)
        {
            if(fabs(z/G.Z_GasDisk[i])>30)
                continue;
        
            double trunc_fac = GetTruncGas(r, i);
            
            if (trunc_fac==0)
                continue;

            double poly_const = GetPolyConst(r, i);
            double dens_const = GetDensConst(r, i);

            //double Sigma_Profile[2], Rho_Profile[2];

            //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

            double rho_0 = dens_const;//*Sigma_Profile[0];

            density += rho_0*exp(psi/poly_const)*trunc_fac;//*Sigma_Profile[0];
        }
        else
        {
            double trunc_fac = GetTruncGas(r, i);

            if (trunc_fac==0)
                continue;

            double poly_const = GetPolyConst(r, i);

            double exponent = 1/(G.Gamma[i]-1);
            double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);

		    double Sigma_Profile[2], Rho_Profile[2];

            GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

            density += pow(-constant*psi, exponent)*Sigma_Profile[0]*trunc_fac;
            //density += pow(-constant*psi, exponent)*trunc_fac;

            //cout << "Gasd " << r << " " << pow(-constant*psi, exponent) << " " 
            //     << "       " << psi << " " << constant << " " << exponent << " " 
            //     << polytrope_const << endl;
        }
    }
        
    return density;   
}

double GasDiskDensf2(double r, double z)
{
    double psi = Pot(r, z);
    return GasDiskDens2(r, z, psi);
}

double GasDiskDensI2(double r, double z, double psi, int i, 
                     double poly_const, double dens_const)
{
    double psi_zero = Pot(r, 0);
    psi -= psi_zero;
    
    double trunc_fac = GetTruncGas(r, i);

    if (trunc_fac==0)
    {
        return 0;
    }
    else if (G.Gamma[i] == 1)
    {
	    //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density, rho_0 = dens_const;//*Sigma_Profile[0];

        density = rho_0*exp(psi/poly_const)*trunc_fac;//*Sigma_Profile[0];

        if (density < 0)
        {
            return 0;
        }

        return density;
    }
    else
    {
        double exponent = 1/(G.Gamma[i]-1);
        double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);

        double Sigma_Profile[2], Rho_Profile[2];

        GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density = pow(-constant*psi, exponent)*Sigma_Profile[0]*trunc_fac;
        //density += pow(-constant*psi, exponent)*trunc_fac;
        
        return density;
    }
}

double GasDiskDensfI2(double r, double z, int i, double poly_const, 
                      double dens_const)
{
    double psi = Pot(r, z);
    return GasDiskDensI2(r, z, psi, i, poly_const, dens_const);
}

double GasDiskDensfI2(double r, double z, int i)
{
    double psi = Pot(r, z);
    double dens_const = GetDensConst(r, i);
    double poly_const = GetPolyConst(r, i);
    return GasDiskDensI2(r, z, psi, i, poly_const, dens_const);
}

double GasDiskDensINoTrunc2(double r, double z, double psi, int i, 
                            double dens_const, double psi_zero)
{
    //double psi_zero = Pot(r, 0);
    psi -= psi_zero;
    
    if (G.Gamma[i] == 1)
    {
            double trunc_fac = GetTruncGas(r, i);
            
            if (trunc_fac==0)
                return 0;

        double poly_const = GetPolyConst(r, i);
        
	    //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density, rho_0 = dens_const;//*Sigma_Profile[0];

        density = rho_0*exp(psi/poly_const)*trunc_fac;
        //pow(-constant*psi, exponent)*trunc_fac;//*Sigma_Profile[0];

        if (density <= 0)
        {
            return 0;
        }

        return density;
    }
//     else
//     {
//         double exponent = 1/(G.Gamma[i]-1);
//         double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);
// 
// 		double Sigma_Profile[2], Rho_Profile[2];
// 
//         GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
// 
//         density += pow(-constant*psi, exponent)*Sigma_Profile[0]*trunc_fac;
//         //density += pow(-constant*psi, exponent)*trunc_fac;
// 
//         //cout << "Gasd " << r << " " << pow(-constant*psi, exponent) << " " 
//         //     << "       " << psi << " " << constant << " " << exponent << " " 
//         //     << polytrope_const << endl;
//     }   
}

double GasDiskDensINoTrunc(double r, double z, double psi, int i, 
                           double poly_const)
{
//     //if (psi < 0)
//     //    return 0;
//    
    if (G.Gamma[i] == 1)
    {
        cout << "Gamma = 1 but the general density function is being called! "
             << "This should not happen! Goodbye." << endl;
        exit(1);
    }
     
    double psi_edge = Pot(r, Radius[nr-1]);
    psi -= psi_edge;
    
    double density;
    
    double exponent = 1/(G.Gamma[i]-1);
    double constant = (1-G.Gamma[i])/(G.Gamma[i]*poly_const);

	double Sigma_Profile[2], Rho_Profile[2];

    GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

    density = pow(-constant*psi, exponent);//*trunc_fac;//*Sigma_Profile[0];
    
    //if (density == 0) 
    //    cout << "dens zero " << r << " " << z << " " << constant << " " << psi 
    //         << "                  " << pow(-constant*psi, exponent) << " " 
    //         << Sigma_Profile[0] << " " << trunc_fac << " " << density << endl;
    
    if (density < 0)
    {
        return 0;
    }
        
    return density;   
}
