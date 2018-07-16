//This file is to calculate the density profile and derivatives of the density
//for the disk. It is designed to be easily modified - comment out the profiles
//you don't need and add in the ones you need. 

#include "galaxy.h"

//function here to calculate disk constants *once*, assuming the diskconsts are 
//global variables defined in galaxy.h

//structure of the file in.diskpars:
// <whitespace delineated parameters> 0
// <whitespace delineated parameters> 0
//etc, where the 0s indicate the end of the parameter list for each disk 
//component. Same structure for the file in.diskdf.

//Define the density and derivative variables first, and initialize to 0. Then 
//calculate the needed densities and derivatives using +=. This makes it easier
//to handle mix and match situations, such as an inner truncated thin disk and
//regular thick disk, for example.

//Edit: This code deals with the gas components. There's no change essentially
//from the stellar disk components since the radial dependence can still be 
//exponential while the vertical dependence (which certainly does change from 
//stellar) is accounted for in the file diskdens.cpp.

void GasDiskProfile(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0;//, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        //ddensity = -density / G.R_GasDisk[i];
        //dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
        
        //double height = 1.8*tanh((radius-38.0)/15)+1.9;
        //rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/height))*cosh(fabs(z/height)));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        //ddensity = -density / G.R_GasDisk[i];
        //dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    //Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
   
void GasDiskProfilePrime(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0;//, dddensity = 0;
    double rho_density = 0;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        ddensity = -density / G.R_GasDisk[i];
        //dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        ddensity = -density / G.R_GasDisk[i];
        //dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
   
void GasDiskProfile2Prime(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        ddensity = -density / G.R_GasDisk[i];
        dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_GasDisk[i]);
        ddensity = -density / G.R_GasDisk[i];
        dddensity = -ddensity / G.R_GasDisk[i];
        rho_density = exp(-radius/G.R_GasDisk[i])/(cosh(fabs(z/G.Z_GasDisk[i]))*cosh(fabs(z/G.Z_GasDisk[i])));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    Sigma_Profile[1] = ddensity;
    Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
   
//in.diskpars should delineate components with a zero. No string libraries
//used here because of the overhead they incur.     
void GetGasDiskParameters(void)
{
    ifstream diskin("in.gasdiskpars", ios::in);
    
    if (!diskin.is_open())
    {
        cerr << "in.gasdiskpars missing. If you don't want a gas disk you need to set\n"
             << "gasdisk_flag to 0 in in.dbh. Exiting..." << endl;
        exit(1);
    }
    
    int param_count = 0;
    double x;
	gas_disk = 0;
    
    cout << "Getting gas disk parameters..." << endl;
    
    while(diskin >> x)
    {
        if (x==-1)
        {
            break;
        }
        else if (x == 0)
        {
            param_count = 0;

            double temp = G.R_GasDisk[gas_disk] * G.R_GasDisk[gas_disk];          
            GasDisk_Const.push_back(G.M_GasDisk[gas_disk] * oneover2pi / temp);
            Rho_GasDisk_Const.push_back(GasDisk_Const[gas_disk] * 0.5 / G.Z_GasDisk[gas_disk]);
            //Polytrope_Gas_Const.push_back(1e2);
            
            ++gas_disk;
            continue;
        }
        else
        {
            ++param_count;
            if (param_count == 1)
                G.M_GasDisk.push_back(x);
            else if (param_count == 2)
                G.R_GasDisk.push_back(x);
            else if (param_count == 3)
                G.Out_GasDisk.push_back(x);
            else if (param_count == 4)
                G.Z_GasDisk.push_back(x);
            else if (param_count == 5)
                G.Dr_Trunc_Gas.push_back(x);
            else if (param_count == 6)
                G.Sigma_0_Gas.push_back(x);
            else if (param_count == 7)
                G.R_Sigma_Gas.push_back(x);
            else if (param_count == 8)
                G.Gamma.push_back(x);
            else if (param_count == 9)
                G.R_Kormendy_Gas.push_back(x);
            else if (param_count == 10)
                G.Alpha_Gas.push_back(x);
            else
                G.Extra_GasDisk_Parameters.push_back(x);
        }
    }
    
    cout << "  Gas Disk: " << G.M_GasDisk[0] << " " << G.R_GasDisk[0] << " " << G.Out_GasDisk[0]
         << " " << G.Z_GasDisk[0] << " " << G.Dr_Trunc_Gas[0] << " " << G.Sigma_0_Gas[0] 
         << " " << G.R_Sigma_Gas[0] << " " << G.Gamma[0] << endl;
}

// void GetDiskDFParameters(void)
// {
//     ifstream diskin("in.diskdf", ios::in);
//     
//     disk_df = 0;
// 	int param_count = 0;
//     double x;
//     
//     while(cin >> x)
//     {
//         if (x == 0)
//         {
//             param_count = 0;
//             ++disk_df;
//             continue;
//         }
//         else
//         {
//             ++param_count;
//             if (param_count == 1)
//                 G.Sigma_0.push_back(x);
//             else if (param_count == 2)
//                 G.R_Sigma.push_back(x);
//             else
//                 G.Extra_DiskDF_Parameters.push_back(x);
//         }
//     }
// }
