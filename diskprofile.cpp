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

void DiskProfile(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0;//, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = G.Z_Disk[i];//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_Disk[i]);
        //ddensity = -density / G.R_Disk[i];
        //dddensity = -ddensity / G.R_Disk[i];
        //rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_Disk[i]);
        //ddensity = -density / G.R_Disk[i];
        //dddensity = -ddensity / G.R_Disk[i];
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    //Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
   
void DiskProfilePrime(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0;//, dddensity = 0;
    double rho_density = 0;
    double zfac = G.Z_Disk[i];//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_Disk[i]);
        ddensity = -density / G.R_Disk[i];
        //dddensity = -ddensity / G.R_Disk[i];
        //rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_Disk[i]);
        ddensity = -density / G.R_Disk[i];
        //dddensity = -ddensity / G.R_Disk[i];
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
   
void DiskProfile2Prime(double &radius, double &z, int &i, double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = G.Z_Disk[i];//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/G.R_Disk[i]);
        ddensity = -density / G.R_Disk[i];
        dddensity = -ddensity / G.R_Disk[i];
        //rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/G.R_Disk[i]);
        ddensity = -density / G.R_Disk[i];
        dddensity = -ddensity / G.R_Disk[i];
        rho_density = exp(-radius/G.R_Disk[i])/(cosh(fabs(z/G.Z_Disk[i]))*cosh(fabs(z/G.Z_Disk[i])));
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
void GetDiskParameters(void)
{
    ifstream diskin("in.diskpars", ios::in);
    
    if (!diskin.is_open())
    {
        cerr << "in.diskpars missing. If you don't want a disk you need to set\n"
             << "disk_flag to 0 in in.dbh. Exiting..." << endl;
        exit(1);
    }
    
    int param_count = 0;
    double x;
	disk = 0;
    
    cout << "Getting disk parameters..." << endl;
    
    while(diskin >> x)
    {
        if (x==-1)
        {
            break;
        }
        else if (x == 0)
        {
            param_count = 0;

            double temp = G.R_Disk[disk] * G.R_Disk[disk];          
            Disk_Const.push_back(G.M_Disk[disk] * oneover2pi / temp);
            Rho_Disk_Const.push_back(Disk_Const[disk] * 0.5 / G.Z_Disk[disk]);
            
            ++disk;
            continue;
        }
        else
        {
            ++param_count;
            if (param_count == 1)
                G.M_Disk.push_back(x);
            else if (param_count == 2)
                G.R_Disk.push_back(x);
            else if (param_count == 3)
                G.Out_Disk.push_back(x);
            else if (param_count == 4)
                G.Z_Disk.push_back(x);
            else if (param_count == 5)
                G.Dr_Trunc.push_back(x);
            else if (param_count == 6)
                G.Sigma_0.push_back(x);
            else if (param_count == 7)
                G.R_Sigma.push_back(x);
            else if (param_count == 8)
                G.R_Kormendy.push_back(x);
            else if (param_count == 9)
                G.Alpha.push_back(x);
            else
                G.Extra_Disk_Parameters.push_back(x);
        }
    }
    
    cout << "  Disk: " << G.M_Disk[0] << " " << G.R_Disk[0] << " " << G.Out_Disk[0]
         << " " << G.Z_Disk[0] << " " << G.Dr_Trunc[0] << " " << G.Sigma_0[0] 
         << " " << G.R_Sigma[0] << endl;
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
