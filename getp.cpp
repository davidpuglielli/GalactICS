#include "galaxy.h"

void CalculatePolytropeConstants(void)
{
    for (int j = 0; j < gas_disk; ++j)
    {
        if (G.Gamma[j] == 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];

                Polytrope_Const[j][i] = G.Sigma_0_Gas[j]*G.Sigma_0_Gas[j];

                if (i%1000==0)
                {
                    double psi = Pot(rad, 0);
                    double rho0 = GasDiskDensINoTrunc2(rad,0,psi,j,Polytrope_Const[j][i], psi);
                    double rho1 = DiskDens(rad,0,psi);
                    double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << Polytrope_Const[j][i] << "               "
                         << psi << " " << rho0 << " " << rho1 << " " << rho2 << endl;
                }
            }
        }
        else if (G.Gamma[j] > 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];

                Polytrope_Const[j][i] = IterateGasConst(rad, j);

                if (i%10000==0)
                {
                    double psi = Pot(rad, 0);
                    double rho0 = GasDiskDensINoTrunc(rad,0,psi,j,Polytrope_Const[j][i]);
                    double rho1 = DiskDens(rad,0,psi);
                    double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << Polytrope_Const[j][i] << "               "
                         << psi << " " << rho0 << " " << rho1 << " " << rho2 << endl;
                }
            }
        }
        else
        {
            cout << "Gamma for gas disk " << j << " must be >= 1. Exiting..."
                 << endl;
            exit(1);
        }
    }
}
    
void CalculateGasDensConstants(void)
{
    for (int j = 0; j < gas_disk; ++j)
    {
        if (G.Gamma[j] == 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];

                GasDensity_Const[j][i] = IterateGasDens(rad, j);

                if (i%1000==0)
                {
                    //cout << "here " << i << j << endl;
                    double psi = Pot(rad, 0);
                    double rho0 = GasDiskDensINoTrunc2(rad,0,psi,j,GasDensity_Const[j][i], psi);
                    double rho1 = DiskDens(rad,0,psi);
                    double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    //cout << "here " << i << j << endl;
                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << GasDensity_Const[j][i] << "               "
                         << Rho_GasDisk_Const[j] << " " 
                         << GasDensity_Const[j][i]/exp(-rad/G.R_GasDisk[j]) 
                         << " " << rho1 << " " << rho2 << endl;
                }
            }
            
            //cout << "Done" << endl;
        }
        else if (G.Gamma[j] > 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];
                double psi = Pot(rad,0);

                GasDensity_Const[j][i] = GasDiskDensINoTrunc(rad,0,psi,j,Polytrope_Const[j][i]);

                if (i%1==0)
                {
                    double psi = Pot(rad, 0);
                    double rho0 = GasDiskDensINoTrunc2(rad,0,psi,j,GasDensity_Const[j][i], psi);
                    double rho1 = DiskDens(rad,0,psi);
                    double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << GasDensity_Const[j][i] << "               "
                         << psi << " " << rho0 << " " << rho1 << " " << rho2 << endl;
                }
            }
        }
        else
        {
            cout << "Gamma for gas disk " << j << " must be >= 1. Exiting..."
                 << endl;
            exit(1);
        }
    }
}
    
void CalculateGasScaleHeights(void)
{
    for (int j = 0; j < gas_disk; ++j)
    {
        if (G.Gamma[j] == 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];//GasDensity_Height[j][i] = G.Z_GasDisk[j];continue;
                double poly = GetDensConst(rad, j);
                double surfden = GasDiskSurfaceDensfI2(rad, j, poly);

                for (double z = 0; z < 10; z+=dr)
                {
                    double psi = Pot(rad,z);
                    double partialsurfden = GasDiskSurfaceDensfI2(rad, j, poly, z);
                    double ratio = partialsurfden/surfden;

                    if(ratio>0.75)
                    {
                        GasDensity_Height[j][i] = max(0.05, z);
                        //cout << "    " << r << " " << setw(12) << z << " "
                        //     << GasDiskDensINoTrunc2(r,z,psi,0,poly) << " "
                        //     << partialsurfden << " " << ratio << "       "
                        //     << DiskDens(r,z,psi) << endl;
                        break;
                    }
                }
        
                if (i%1000==0)
                {
                    //double psi = Pot(rad, 0);
                    //double rho0 = GasDiskDensINoTrunc2(rad,0,psi,j,GasDensity_Const[j][i]);
                    //double rho1 = DiskDens(rad,0,psi);
                    //double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << GasDensity_Height[j][i] << endl;
                }
            }
        }
        else if (G.Gamma[j] > 1)
        {
            for (int i = 0; i < nr; ++i)
            {
                double rad = Radius[i];
                double psi = Pot(rad,0);

                GasDensity_Const[j][i] = GasDiskDensINoTrunc(rad,0,psi,j,Polytrope_Const[j][i]);

                if (i%10000==0)
                {
                    double psi = Pot(rad, 0);
                    double rho0 = GasDiskDensINoTrunc2(rad,0,psi,j,GasDensity_Const[j][i], psi);
                    double rho1 = DiskDens(rad,0,psi);
                    double rho2 = GasDisk_Const[j]*exp(-rad/G.R_GasDisk[j]);

                    cout << "Did " << i << " calculation " << setw(12) << rad 
                         << " " << GasDensity_Const[j][i] << "               "
                         << psi << " " << rho0 << " " << rho1 << " " << rho2 << endl;
                }
            }
        }
        else
        {
            cout << "Gamma for gas disk " << j << " must be >= 1. Exiting..."
                 << endl;
            exit(1);
        }
    }
}
    
double GetPolyConst(double r, int j)
{
    if (G.Gamma[j] == 1)
    {
        return Polytrope_Const[j].at(0);
    }
    else
    {
        int ihi = ceil(r/dr);

        if(r < dr)
        {
            ihi = 1;
        }
        else if (ihi < 1)
        {
            cout << "GetPolyConst finds out of range indices. Exiting..." << endl;
            exit(1);
        }
        else if (ihi > nr-1)
        {
            ihi = nr-1;
        }

        double r1 = Radius[ihi-1];
        double r2 = Radius[ihi];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << r << " " << ihi << " " << t << " " << Polytrope_Const[j][ihi] << " "
        //     << tm1 << " " << Polytrope_Const[j][ihi-1] << endl;
        return t*Polytrope_Const[j][ihi]+tm1*Polytrope_Const[j][ihi-1];
    }
}    

double GetDensConst(double r, int j)
{
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetDensConst finds out of range indices. Exiting..." << endl;
        exit(1);
    }
    else if (ihi > nr-1)
    {
        ihi = nr-1;
    }
    
    double r1 = Radius[ihi-1];
    double r2 = Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    //cout << r << " " << ihi << " " << t << " " << Polytrope_Const[j][ihi] << " "
    //     << tm1 << " " << Polytrope_Const[j][ihi-1] << endl;
    return t*GasDensity_Const[j][ihi]+tm1*GasDensity_Const[j][ihi-1];
}    
    
double GetGasHeight(double r, int j)
{
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetDensConst finds out of range indices. Exiting..." << endl;
        exit(1);
    }
    else if (ihi > nr-1)
    {
        ihi = nr-1;
    }
    
    double r1 = Radius[ihi-1];
    double r2 = Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    //cout << r << " " << ihi << " " << t << " " << Polytrope_Const[j][ihi] << " "
    //     << tm1 << " " << Polytrope_Const[j][ihi-1] << endl;
    return t*GasDensity_Height[j][ihi]+tm1*GasDensity_Height[j][ihi-1];
}    
    
double GasDiskSurfaceDensfI(double r, int j, double poly_const)
{
    double dz = 0.01;
    double zmax = 5;
    int integ_pts = int(zmax/dz);  
    double psi = Pot(r, 0);
    
    double plane_density = GasDiskDensINoTrunc(r, 0, psi, j, poly_const)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    
    for (int i = 1; i < integ_pts; ++i)
    {
        psi = Pot(r, i*dz);
        //if(poly_const<10)
        //cout << "iterate1 " << r << " " << i*dz << " " << density << " " 
        //     << psi << endl;
        density += GasDiskDensINoTrunc(r, i*dz, psi, j, poly_const)*dz;
    }
    
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << poly_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
double GasDiskSurfaceDensfI2(double r, int j, double dens_const)
{
    double dz = 0.01;
    double zmax = 5;
    int integ_pts = int(zmax/dz);  
    double psi = Pot(r, 0);
    double psi_zero = psi;
    
    double plane_density = GasDiskDensINoTrunc2(r, 0, psi, j, dens_const, psi_zero)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    
    for (int i = 1; i < integ_pts; ++i)
    {
        psi = Pot(r, i*dz);
        //if(poly_const<10)
        //cout << "iterate1 " << r << " " << i*dz << " " << density << " " 
        //     << psi << endl;
        density += GasDiskDensINoTrunc2(r, i*dz, psi, j, dens_const, psi_zero)*dz;
    }
    
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << dens_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
double GasDiskSurfaceDensfI(double r, int j, double poly_const, double zmax)
{
    double dz = 0.01;
    int integ_pts = int(zmax/dz);  
    double psi = Pot(r, 0);
    
    double plane_density = GasDiskDensINoTrunc(r, 0, psi, j, poly_const)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    
    for (int i = 1; i < integ_pts; ++i)
    {
        psi = Pot(r, i*dz);
        //if(poly_const<10)
        //  cout << "iterate1 " << r << " " << i*dz << " " 
        //       << density << " " << psi << endl;
        density += GasDiskDensINoTrunc(r, i*dz, psi, j, poly_const)*dz;
    }
    
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << poly_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
double GasDiskSurfaceDensfI2(double r, int j, double dens_const, double zmax)
{
    double dz = 0.01;
    int integ_pts = int(zmax/dz);  
    double psi = Pot(r, 0);
    double psi_zero = psi;
    
    double plane_density = GasDiskDensINoTrunc2(r, 0, psi, j, dens_const, psi_zero)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    
    for (int i = 1; i < integ_pts; ++i)
    {
        psi = Pot(r, i*dz);
        //if(poly_const<10)
        //cout << "iterate1 " << r << " " << i*dz << " " << density << " " 
        //     << psi << endl;
        density += GasDiskDensINoTrunc2(r, i*dz, psi, j, dens_const, psi_zero)*dz;
    }
    
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << poly_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
double IterateGasConst(double r, int j)
{
    double Sigma_Profile[2], Rho_Profile[2], zz = 0;
		
    GasDiskProfile(r, zz, j, Sigma_Profile, Rho_Profile);

    double target_density = GasDisk_Const[j]*Sigma_Profile[0]; 
       
    //Find brackets for the polytropic constant. Here K is a function of density
    //so the following algorithm is a transcription of the root bisection
    //found in gendf.cpp. First, the find brackets part
    double polyconst_lower, polyconst_upper;
    polyconst_lower = polyconst_upper = 100;
    
    double rho_upper = GasDiskSurfaceDensfI(r, j, polyconst_upper);
    double rho_lower = GasDiskSurfaceDensfI(r, j, polyconst_lower);
    
    //cout << "                                    here" << endl;
    while (rho_lower <= target_density)
    {
        polyconst_lower *= 10;//0.1;
        rho_lower = GasDiskSurfaceDensfI(r, j, polyconst_lower);
        //if (polyconst_lower<100)
            //cout << "loop " << r << " " << polyconst_lower << " " 
            //     << rho_lower << " " << target_density <<  endl;
        
        if (rho_lower == 0)
            return 0;
    }
    
    //cout << "                                    here" << endl;
    while (rho_upper >= target_density)
    {
        polyconst_upper *= 0.1;//10;
        rho_upper = GasDiskSurfaceDensfI(r, j, polyconst_upper);
        //cout << r << " " << polyconst_upper << " " << rho_upper << " " 
        //     << target_density << endl;
    }
    
    //cout << "                                    here2" << endl;
    
    //Now the root bisection part
    rho_upper = GasDiskSurfaceDensfI(r, j, polyconst_upper) - target_density;
    rho_lower = GasDiskSurfaceDensfI(r, j, polyconst_lower) - target_density;
    
    double polyconst_midpoint, rho_midpoint, poly_gap;
    int iters = 0;
    
    if (rho_upper*rho_lower > 0)
    {
        cerr << "Bisection endpoints do not bracket root. Exiting..." << endl;
        cerr << rho_upper << " " << rho_lower << " " << polyconst_lower << " " 
             << polyconst_upper << endl;
        exit(1);
    }
    
    do
    {
        ++iters;
        
        polyconst_midpoint = (polyconst_upper+polyconst_lower)*0.5;
        
        rho_midpoint = GasDiskSurfaceDensfI(r, j, polyconst_midpoint) - 
                       target_density;
        
        if (rho_lower*rho_midpoint < 0)
        {
            polyconst_upper = polyconst_midpoint;
            rho_upper = rho_midpoint;
        }
        else if (rho_upper*rho_midpoint < 0)
        {
            polyconst_lower = polyconst_midpoint;
            rho_lower = rho_midpoint;
        }
        else if (rho_midpoint==0)
        {
            return polyconst_midpoint;
            cout << "\nFinished root bisection\n" << endl;
        }
        else
        {
            cerr << "Root bisection for polyconst failed! exiting..." << endl;
            cerr << rho_lower << " " << rho_upper << " " << rho_midpoint << " "
                 << polyconst_lower << " " << polyconst_upper << " " 
                 << polyconst_midpoint << " " << target_density << endl;
            exit(1);
        }
        
        poly_gap = polyconst_upper-polyconst_lower;
    }
    while (iters < 200 && poly_gap > 1e-6*polyconst_midpoint);
    
    //cout << "                                    " << iters << " " 
    //     << target_density << " " << rho_midpoint << endl;
    
    return polyconst_midpoint;
}
        
double IterateGasDens(double r, int j)
{
    //cout << "iterate " << endl;
    
    double Sigma_Profile[2], Rho_Profile[2], zz = 0;
		
    GasDiskProfile(r, zz, j, Sigma_Profile, Rho_Profile);

    double target_density = GasDisk_Const[j]*Sigma_Profile[0]; 
       
    //cout << "iterate " << endl;
    //Find brackets for the polytropic constant. Here K is a function of density
    //so the following algorithm is a transcription of the root bisection
    //found in gendf.cpp. First, the find brackets part
    double dens_lower, dens_upper;
    dens_lower = dens_upper = 0.1;
    
    //cout << "iterate " << endl;
    double rho_upper = GasDiskSurfaceDensfI2(r, j, dens_upper);
    //cout << "iterate " << endl;
    double rho_lower = GasDiskSurfaceDensfI2(r, j, dens_lower);
    
    //cout << "iterate " << endl;
    while (rho_lower >= target_density)
    {
        dens_lower *= 0.1;
        rho_lower = GasDiskSurfaceDensfI2(r, j, dens_lower);
        
        if (rho_lower == 0)
            return 0;
    }
    
    while (rho_upper <= target_density)
    {
        dens_upper *= 10;
        rho_upper = GasDiskSurfaceDensfI2(r, j, dens_upper);
        
        if (rho_upper == 0)
            return 0;
    }
    
    //Now the root bisection part
    rho_upper = GasDiskSurfaceDensfI2(r, j, dens_upper) - target_density;
    rho_lower = GasDiskSurfaceDensfI2(r, j, dens_lower) - target_density;
    
    double dens_midpoint, rho_midpoint, poly_gap;
    int iters = 0;
    
    if (rho_upper*rho_lower > 0)
    {
        cerr << "Bisection endpoints do not bracket root. Exiting..." << endl;
        cerr << rho_upper << " " << rho_lower << " " << dens_lower << " " 
             << dens_upper << endl;
        exit(1);
    }
    
    do
    {
        ++iters;//cout<<"here "<<endl;
        
        dens_midpoint = (dens_upper+dens_lower)*0.5;
        
        rho_midpoint = GasDiskSurfaceDensfI2(r, j, dens_midpoint) - 
                       target_density;
        
        if (rho_lower*rho_midpoint < 0)
        {
            dens_upper = dens_midpoint;
            rho_upper = rho_midpoint;
        }
        else if (rho_upper*rho_midpoint < 0)
        {
            dens_lower = dens_midpoint;
            rho_lower = rho_midpoint;
        }
        else if (rho_midpoint==0)
        {
            return dens_midpoint;
            cout << "\nFinished root bisection\n" << endl;
        }
        else
        {
            cerr << "Root bisection for dens failed! exiting..." << endl;
            cerr << rho_lower << " " << rho_upper << " " << rho_midpoint << " "
                 << dens_lower << " " << dens_upper << " " 
                 << dens_midpoint << " " << target_density << endl;
            exit(1);
        }
        
        poly_gap = dens_upper-dens_lower;
    }
    while (iters < 200 && poly_gap > 1e-2*dens_midpoint);
    
    //cout << "                                    " << iters << " " 
    //     << target_density << " " << dens_midpoint << endl;
    
    return dens_midpoint;
}
        
// void GetP(void)
// {
//     for (int i = 0; i < Radius.size(); ++i)
//     {
//         double r = Radius[i];
//         
//         double psi_zh = Pot(r, G.Z_Disk[0]);
//         double psi_0 = Pot(r, 0.0);
//         double true_stellar_sig_z2 = (psi_zh-psi_0)/log(0.419974);
//         double gas_sig_z2 = G.Z_GasDisk[0]
// 
//         double p = true_stellar_sig_z2/gas_sig_z2;
//         
//         P_Ratio[i] = p;
//     }
// }
