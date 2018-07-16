//Approximate disk density, force and potential
//Note I still have to sort out r, s and radius variables

#include "galaxy.h"

void AppDiskForce(double &s, double &z, double &fsad, double &fzad)
{
    fsad=0;fzad=0;
    
    double r = sqrt(s*s + z*z);
    fsad = fzad = 0;
        
    for (int i = 0; i < disk; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrime(r, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[2], Rho_Profile[2];

        DiskProfilePrime(r, z, i, Sigma_Profile, Rho_Profile);
    
        double zz = fabs(z/G.Z_Disk[i]);
        double e2zz = exp(-2*zz);
        double tlncoshz = zz+log(0.5*(1.0 + e2zz));
        
		double f = Sigma_Profile[0]*trunc_fac;
		double f_prime = Sigma_Profile[0]*trunc_facprime + 
                         Sigma_Profile[1]*trunc_fac;
		
        fsad += -fourpi*Rho_Disk_Const[i]*f_prime*s*tlncoshz;
        fzad += -fourpi*Rho_Disk_Const[i]*
                (f_prime*z*tlncoshz + f/G.Z_Disk[i]*tanh(z/G.Z_Disk[i]));
    }

//     double r = s;//sqrt(s*s+z*z);
//     fsad = 0; fzad = 0;
//     
//     if (r == 0 && z == 0)
//     {
//         for (int i = 0; i < disk; ++i)
//         {
//             fsad = App_Disk_Force_R[i][0][0];
//             fzad = App_Disk_Force_Z[i][0][0];
//             return;
//         }
//     }
//     
//     int iri = ceil(r/dr);
//     int izi = ceil(fabs(z)/dr);
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     for (int i = 0; i < disk; ++i)
//     {
//         double force1 = t*App_Disk_Force_R[i][iri][izi] + tm1*App_Disk_Force_R[i][iri-1][izi];
//         double force2 = t*App_Disk_Force_R[i][iri][izi-1] + tm1*App_Disk_Force_R[i][iri-1][izi-1];
//         double force3 = ss*App_Disk_Force_R[i][iri][izi] + sm1*App_Disk_Force_R[i][iri][izi-1];
//         double force4 = ss*App_Disk_Force_R[i][iri-1][izi] + sm1*App_Disk_Force_R[i][iri-1][izi-1];
//     
//         fsad += 0.25*(force1+force2+force3+force4);
//         
//         force1 = t*App_Disk_Force_Z[i][iri][izi] + tm1*App_Disk_Force_Z[i][iri-1][izi];
//         force2 = t*App_Disk_Force_Z[i][iri][izi-1] + tm1*App_Disk_Force_Z[i][iri-1][izi-1];
//         force3 = ss*App_Disk_Force_Z[i][iri][izi] + sm1*App_Disk_Force_Z[i][iri][izi-1];
//         force4 = ss*App_Disk_Force_Z[i][iri-1][izi] + sm1*App_Disk_Force_Z[i][iri-1][izi-1];
//     
//         fzad += 0.25*(force1+force2+force3+force4);
//     }
}
        
double AppDiskPot(double &s, double &z)
{
    //return 0;
    
    double appdiskpot = 0;
    double radius = sqrt(s*s+z*z);
    
    //Now add up the appdiskpot contributions from each sech^2 component
    //profile[3*1] is the density of the i-th component
    for (int i = 0; i < disk; ++i)
    {
        double Sigma_Profile[2], Rho_Profile[2];

        DiskProfile(radius, z, i, Sigma_Profile, Rho_Profile);
        
        double trunc_fac = GetTrunc(radius, i);
        double zz = fabs(z/G.Z_Disk[i]);
        appdiskpot += -fourpi*Rho_Disk_Const[i]*Sigma_Profile[0]*G.Z_Disk[i]*
                      G.Z_Disk[i]*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
    }
    
    return appdiskpot;

//     double r = s;//sqrt(s*s+z*z);
//     double appdiskpot = 0;
//     
//     //if (r == 0 && z == 0)
//     //{
//     //    return App_GasDisk_Dens;
//     //}
//     
//     int iri = ceil(r/dr);    
//     int izi = ceil(fabs(z)/dr);    
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     for (int i = 0; i < disk; ++i)
//     {
//         double pot1 = t*App_Disk_Pot[i][iri][izi] + tm1*App_Disk_Pot[i][iri-1][izi];
//         double pot2 = t*App_Disk_Pot[i][iri][izi-1] + tm1*App_Disk_Pot[i][iri-1][izi-1];
//         double pot3 = ss*App_Disk_Pot[i][iri][izi] + sm1*App_Disk_Pot[i][iri][izi-1];
//         double pot4 = ss*App_Disk_Pot[i][iri-1][izi] + sm1*App_Disk_Pot[i][iri-1][izi-1];
//     
//         appdiskpot += ss*pot1+sm1*pot2;
//         //appdiskpot += 0.25*(pot1+pot2+pot3+pot4);
//     }
//    
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     return appdiskpot;
}

double AppDiskDens(double &s, double &z)
{
    //return 0;
    
// KK: This is the density corresponding to the first-guess disk potential
// f(r)*erfc((r-outdisk)/sqrt(2)drtrunc)/2 * 4 pi G zdisk**2 log(z/zdisk)
// where r is spherical radius. f(r) is here taken as an exponential.
//
// The corresponding density is:
// f(r)*erfc*sech(z/zdisk)**2 + radial gradient terms.
// For radii below one scale radius, we have replaced the radial exponential 
// (which gives rise to a singular laplacian) with a quartic that joins on 
// smoothly.

// DP: The last part of KK's comment appears to never have been implemented.

// DP: See equation 12 of Kuijken & Dubinski (1995) for the return value
    
    double r = sqrt(s*s+z*z);
    double appdiskdens = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrime(r, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[3], Rho_Profile[3];

        DiskProfile2Prime(r, z, i, Sigma_Profile, Rho_Profile);
        
        double f = Sigma_Profile[0]*G.Z_Disk[i]*G.Z_Disk[i]*trunc_fac;
        double f1r = G.Z_Disk[i]*G.Z_Disk[i]*
                     (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/r;
        double f2 = G.Z_Disk[i]*G.Z_Disk[i]*(Sigma_Profile[2]*trunc_fac +
                    2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		            (r-G.Out_Disk[i])/G.Dr_Trunc[i]/G.Dr_Trunc[i]);
        
        if (r==0)
        {
            f1r = f2 = 0;
        }        
    
        double zz = fabs(z/G.Z_Disk[i]);
        double ezz = exp(-zz), e2zz = exp(-2*zz);
        double tlncosh = zz+log(0.5*(1.0 + e2zz));
        double tztanh = zz*(1-e2zz)/(1+e2zz);
        double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
        double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
                       f*tsech2/G.Z_Disk[i]/G.Z_Disk[i];
        appdiskdens += Rho_Disk_Const[i]*total;
        //cout<<"appdiskdens "<<Rho_Disk_Const[i]*f<<" "<<Rho_Disk_Const[i]*f1r
        //    <<" "<<Rho_Disk_Const[i]*f2<<" "<<tlncosh<<" "<<tztanh<<" "<<tsech2<<endl;
        //cout<<"appdiskdens1 "<<" "<<f2<<" "<<f1r<<" "<<f<<endl;
        //cout<<"appdiskdens2 "<<" "<<appdiskdens<<endl;
    }
    
    //cout << "dens " << r << " " << z << " " << appdiskdens << endl;
    return appdiskdens;
    
//     double r = s;//sqrt(s*s+z*z);
//     double appdiskdens = 0;
//     
//     //if (r == 0 && z == 0)
//     //{
//     //    return App_Disk_Dens;
//     //}
//     
//     int iri = ceil(r/dr);    
//     int izi = ceil(fabs(z)/dr);    
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     for (int i = 0; i < disk; ++i)
//     {
//         double dens1 = t*App_Disk_Dens[i][iri][izi] + tm1*App_Disk_Dens[i][iri-1][izi];
//         double dens2 = t*App_Disk_Dens[i][iri][izi-1] + tm1*App_Disk_Dens[i][iri-1][izi-1];
//         double dens3 = ss*App_Disk_Dens[i][iri][izi] + sm1*App_Disk_Dens[i][iri][izi-1];
//         double dens4 = ss*App_Disk_Dens[i][iri-1][izi] + sm1*App_Disk_Dens[i][iri-1][izi-1];
//     
//         appdiskdens += ss*dens1+sm1*dens2;
//         //appdiskdens += 0.25*(dens1+dens2+dens3+dens4);
//     }
//     
//     //cout << s << " " << z << " " << appdiskdens << " " << AppDiskDens3(s,z) << endl;
//     
//     return appdiskdens;
}

double AppDiskPot3(double &s, double &z)
{
    //return 0;
    
    double appdiskpot = 0;
    double radius = sqrt(s*s+z*z);
    
    //Now add up the appdiskpot contributions from each sech^2 component
    //profile[3*1] is the density of the i-th component
    for (int i = 0; i < disk; ++i)
    {
        double Sigma_Profile[2], Rho_Profile[2];

        DiskProfile(radius, z, i, Sigma_Profile, Rho_Profile);
        
        double trunc_fac = GetTrunc(radius, i);
        double zz = fabs(z/G.Z_Disk[i]);
        appdiskpot += -fourpi*Rho_Disk_Const[i]*Sigma_Profile[0]*G.Z_Disk[i]*
                      G.Z_Disk[i]*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
    }
    
    return appdiskpot;
}

double AppDiskDens3(double &s, double &z)
{
    //return 0;
    
// KK: This is the density corresponding to the first-guess disk potential
// f(r)*erfc((r-outdisk)/sqrt(2)drtrunc)/2 * 4 pi G zdisk**2 log(z/zdisk)
// where r is spherical radius. f(r) is here taken as an exponential.
// 
// The corresponding density is:
// f(r)*erfc*sech(z/zdisk)**2 + radial gradient terms.
// For radii below one scale radius, we have replaced the radial exponential 
// (which gives rise to a singular laplacian) with a quartic that joins on 
// smoothly.
// 
// DP: The last part of KK's comment appears to never have been implemented.
// 
// DP: See equation 12 of Kuijken & Dubinski (1995) for the return value
    
    double r = sqrt(s*s+z*z);
    double appdiskdens = 0;
    
    for (int i = 0; i < disk; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrime(r, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[3], Rho_Profile[3];

        DiskProfile2Prime(r, z, i, Sigma_Profile, Rho_Profile);
        
        double f = Sigma_Profile[0]*G.Z_Disk[i]*G.Z_Disk[i]*trunc_fac;
        double f1r = G.Z_Disk[i]*G.Z_Disk[i]*
                     (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/r;
        double f2 = G.Z_Disk[i]*G.Z_Disk[i]*(Sigma_Profile[2]*trunc_fac +
                    2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		            (r-G.Out_Disk[i])/G.Dr_Trunc[i]/G.Dr_Trunc[i]);
        
        if (r==0)
        {
            f1r = f2 = 0;
        }        
    
        double zz = fabs(z/G.Z_Disk[i]);
        double ezz = exp(-zz), e2zz = exp(-2*zz);
        double tlncosh = zz+log(0.5*(1.0 + e2zz));
        double tztanh = zz*(1-e2zz)/(1+e2zz);
        double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
        double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
                       f*tsech2/G.Z_Disk[i]/G.Z_Disk[i];
        appdiskdens += Rho_Disk_Const[i]*total;
        //cout<<"appdiskdens "<<Rho_Disk_Const[i]*f<<" "<<Rho_Disk_Const[i]*f1r
        //    <<" "<<Rho_Disk_Const[i]*f2<<" "<<tlncosh<<" "<<tztanh<<" "<<tsech2<<endl;
        //cout<<"appdiskdens1 "<<" "<<f2<<" "<<f1r<<" "<<f<<endl;
        //cout<<"appdiskdens2 "<<" "<<appdiskdens<<endl;
    }
    
    //cout << "dens " << r << " " << z << " " << appdiskdens << endl;
    return appdiskdens;
}

////////////////////////////////////////////////////////////////////////////////    
//Now for the gas disk

void AppGasDiskForce(double &s, double &z, double &fsad, double &fzad)
{
    fsad=0;fzad=0;return;
    
    double r = sqrt(s*s + z*z);
    fsad = fzad = 0;
        
    for (int i = 0; i < gas_disk; ++i)
    {
        double trunc_fac, trunc_facprime;
        double height = GetGasHeight(s, i);
        //double dens_const = GetDensConst(r, i);
        double rho_const = Rho_GasDisk_Const[i];///height;

        GetTruncPrimeGas(r, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[2], Rho_Profile[2];

        GasDiskProfilePrime(s, z, i, Sigma_Profile, Rho_Profile);
    
        double zz = fabs(z/height);
        double e2zz = exp(-2*zz);
        double tlncoshz = zz+log(0.5*(1.0 + e2zz));
        
		double f = Sigma_Profile[0]*trunc_fac;
		double f_prime = Sigma_Profile[0]*trunc_facprime + 
                         Sigma_Profile[1]*trunc_fac;
		
        fsad += -fourpi*rho_const*f_prime*s*tlncoshz;
        fzad += -fourpi*rho_const*
                (f_prime*z*tlncoshz + f/height*tanh(z/height));
    }

//     double r = sqrt(s*s+z*z);
//     fsad = 0; fzad = 0;
//     
//     if (r == 0 && z == 0)
//     {
//         for (int i = 0; i < gas_disk; ++i)
//         {
//             fsad = App_GasDisk_Force_R[i][0][0];
//             fzad = App_GasDisk_Force_Z[i][0][0];
//             return;
//         }
//     }
//     
//     int iri = ceil(r/dr);
//     int izi = ceil(fabs(z)/dr);
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double force1 = t*App_GasDisk_Force_R[i][iri][izi] + tm1*App_GasDisk_Force_R[i][iri-1][izi];
//         double force2 = t*App_GasDisk_Force_R[i][iri][izi-1] + tm1*App_GasDisk_Force_R[i][iri-1][izi-1];
//         double force3 = ss*App_GasDisk_Force_R[i][iri][izi] + sm1*App_GasDisk_Force_R[i][iri][izi-1];
//         double force4 = ss*App_GasDisk_Force_R[i][iri-1][izi] + sm1*App_GasDisk_Force_R[i][iri-1][izi-1];
//     
//         fsad += 0.25*(force1+force2+force3+force4);
//         
//         force1 = t*App_GasDisk_Force_Z[i][iri][izi] + tm1*App_GasDisk_Force_Z[i][iri-1][izi];
//         force2 = t*App_GasDisk_Force_Z[i][iri][izi-1] + tm1*App_GasDisk_Force_Z[i][iri-1][izi-1];
//         force3 = ss*App_GasDisk_Force_Z[i][iri][izi] + sm1*App_GasDisk_Force_Z[i][iri][izi-1];
//         force4 = ss*App_GasDisk_Force_Z[i][iri-1][izi] + sm1*App_GasDisk_Force_Z[i][iri-1][izi-1];
//     
//         fzad += 0.25*(force1+force2+force3+force4);
//     }
}
        
double AppGasDiskPot(double &s, double &z)
{
    return 0;
    
//     double appdiskpot = 0;
//     double radius = sqrt(r*r+z*z);
//     
//     //Now add up the appdiskpot contributions from each sech^2 component
//     //profile[3*1] is the density of the i-th component
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double Sigma_Profile[2], Rho_Profile[2];
// 
//         GasDiskProfile(radius, z, i, Sigma_Profile, Rho_Profile);
//         
//         double trunc_fac = GetTruncGas(radius, i);
//         double zz = fabs(z/G.Z_GasDisk[i]);
//         appdiskpot += -fourpi*Rho_GasDisk_Const[i]*Sigma_Profile[0]*G.Z_GasDisk[i]*
//                       G.Z_GasDisk[i]*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
//     }
//     
//     //cout << "pot " << r << " " << z << " " << appdiskpot << endl;
//     return appdiskpot;

    double appdiskpot = 0;
    double r = sqrt(s*s+z*z);
    
    //Now add up the appdiskpot contributions from each sech^2 component
    //profile[3*1] is the density of the i-th component
    for (int i = 0; i < gas_disk; ++i)
    {
        double height = GetGasHeight(r, i);
        //double dens_const = GetDensConst(r, i);
        double rho_const = Rho_GasDisk_Const[i];///height;
        
        //Let's see if this works
        //double vert_trunc = GetTruncGas(z, i);
        
        double Sigma_Profile[2], Rho_Profile[2];

        GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
        
        double trunc_fac = GetTruncGas(r, i);
        
        if (trunc_fac == 0)
        {
            continue;
        }
        
        double zz = fabs(z/height);
        appdiskpot += -fourpi*rho_const*Sigma_Profile[0]*height*
                      height*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;// * vert_trunc;
    }
    
    //cout << "pot " << r << " " << z << " " << appdiskpot << endl;
    return appdiskpot;

//     double r = sqrt(s*s+z*z);
//     double appdiskpot = 0;
//     
//     //if (r == 0 && z == 0)
//     //{
//     //    return App_GasDisk_Dens;
//     //}
//     
//     int iri = ceil(r/dr);    
//     int izi = ceil(fabs(z)/dr);    
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double pot1 = t*App_GasDisk_Pot[i][iri][izi] + tm1*App_GasDisk_Pot[i][iri-1][izi];
//         double pot2 = t*App_GasDisk_Pot[i][iri][izi-1] + tm1*App_GasDisk_Pot[i][iri-1][izi-1];
//         double pot3 = ss*App_GasDisk_Pot[i][iri][izi] + sm1*App_GasDisk_Pot[i][iri][izi-1];
//         double pot4 = ss*App_GasDisk_Pot[i][iri-1][izi] + sm1*App_GasDisk_Pot[i][iri-1][izi-1];
//     
//         appdiskpot += ss*pot1+sm1*pot2;
//     }
//     
//     //cout << "app " << s << " " << appdiskpot << " " << iri << " " << izi << endl;
//     return appdiskpot;
}

double AppGasDiskDens(double &s, double &z)
{
    return 0;
    
// KK: This is the density corresponding to the first-guess disk potential
// f(r)*erfc((r-outdisk)/sqrt(2)drtrunc)/2 * 4 pi G zdisk**2 log(z/zdisk)
// where r is spherical radius. f(r) is here taken as an exponential.
//
// The corresponding density is:
// f(r)*erfc*sech(z/zdisk)**2 + radial gradient terms.
// For radii below one scale radius, we have replaced the radial exponential 
// (which gives rise to a singular laplacian) with a quartic that joins on 
// smoothly.

// DP: The last part of KK's comment appears to never have been implemented.

// DP: See equation 12 of Kuijken & Dubinski (1995) for the return value
    
//     double r = sqrt(s*s+z*z);
//     double appdiskdens = 0;
//     
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double trunc_fac, trunc_facprime;
// 
//         GetTruncPrimeGas(r, i, trunc_fac, trunc_facprime);
// 
//         double Sigma_Profile[3], Rho_Profile[3];
// 
//         GasDiskProfile2Prime(r, z, i, Sigma_Profile, Rho_Profile);
//         
//         double f = Sigma_Profile[0]*G.Z_GasDisk[i]*G.Z_GasDisk[i]*trunc_fac;
//         double f1r = G.Z_GasDisk[i]*G.Z_GasDisk[i]*
//                      (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/r;
//         double f2 = G.Z_GasDisk[i]*G.Z_GasDisk[i]*(Sigma_Profile[2]*trunc_fac +
//                     2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
// 		            (r-G.Out_GasDisk[i])/G.Dr_Trunc_Gas[i]/G.Dr_Trunc_Gas[i]);
//         
//         if (r==0)
//         {
//             f1r = f2 = 0;
//         }        
//     
//         double zz = fabs(z/G.Z_GasDisk[i]);
//         double ezz = exp(-zz), e2zz = exp(-2*zz);
//         double tlncosh = zz+log(0.5*(1.0 + e2zz));
//         double tztanh = zz*(1-e2zz)/(1+e2zz);
//         double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
//         double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
//                        f*tsech2/G.Z_GasDisk[i]/G.Z_GasDisk[i];
//         appdiskdens += Rho_GasDisk_Const[i]*total;
//         //cout<<"appdiskdens "<<Rho_Disk_Const[i]*f<<" "<<Rho_Disk_Const[i]*f1r
//         //    <<" "<<Rho_Disk_Const[i]*f2<<" "<<tlncosh<<" "<<tztanh<<" "<<tsech2<<endl;
//         //cout<<"appdiskdens1 "<<" "<<f2<<" "<<f1r<<" "<<f<<endl;
//         //cout<<"appdiskdens2 "<<" "<<appdiskdens<<endl;
//     }
//     
//     //cout << "dens " << r << " " << z << " " << appdiskdens << endl;
//     return appdiskdens;

    double r = sqrt(s*s+z*z);
    double appdiskdens = 0;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        double height = GetGasHeight(r, i);
        //double dens_const = GetDensConst(r, i);
        double rho_const = Rho_GasDisk_Const[i];///height;
        
        double trunc_fac, trunc_facprime;

        GetTruncPrimeGas(r, i, trunc_fac, trunc_facprime);

        if (trunc_fac == 0)
        {
            continue;
        }
        
        double Sigma_Profile[3], Rho_Profile[3];

        GasDiskProfile2Prime(r, z, i, Sigma_Profile, Rho_Profile);
        //dens_const /= Sigma_Profile[0];
        
        double f = Sigma_Profile[0]*height*height*trunc_fac;
        double f1r = height*height*
                     (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/r;
        double f2 = height*height*(Sigma_Profile[2]*trunc_fac +
                    2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		            (r-G.Out_GasDisk[i])/G.Dr_Trunc_Gas[i]/G.Dr_Trunc_Gas[i]);
        
        if (r==0)
        {
            f1r = f2 = 0;
        }        
    
        double zz = fabs(z/height);
        double ezz = exp(-zz), e2zz = exp(-2*zz);
        double tlncosh = zz+log(0.5*(1.0 + e2zz));
        double tztanh = zz*(1-e2zz)/(1+e2zz);
        double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
        double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
                       f*tsech2/height/height;
        appdiskdens += rho_const*total;
        //cout<<"appdiskdens "<<Rho_Disk_Const[i]*f<<" "<<Rho_Disk_Const[i]*f1r
        //    <<" "<<Rho_Disk_Const[i]*f2<<" "<<tlncosh<<" "<<tztanh<<" "<<tsech2<<endl;
        //cout<<"appdiskdens1 "<<" "<<f2<<" "<<f1r<<" "<<f<<endl;
        //cout<<"appdiskdens2 "<<" "<<appdiskdens<<endl;
    }
    
    //cout << "dens " << r << " " << z << " " << appdiskdens << endl;
    return appdiskdens;
    
//     double r = sqrt(s*s+z*z);
//     double appdiskdens = 0;
//     
//     //if (r == 0 && z == 0)
//     //{
//     //    return App_GasDisk_Dens;
//     //}
//     
//     int iri = ceil(r/dr);    
//     int izi = ceil(fabs(z)/dr);    
//     
//     if (r < dr)
//     {
//         iri = 1;
//     }
//     else if (iri > nr-1)
//     {
//         iri = nr-1;
//     }
//     
//     
//     if (fabs(z) < dr)
//     {
//         izi = 1;
//     }
//     else if (izi > nr-1)
//     {
//         izi = nr-1;
//     }
//     
//     double r1 = Radius[iri-1];
//     double r2 = Radius[iri];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     double z1 = Radius[izi-1];
//     double z2 = Radius[izi];
//     double ss = (z-z1)/(z2-z1);
//     double sm1 = 1-ss;
//     
//     for (int i = 0; i < gas_disk; ++i)
//     {
//         double dens1 = t*App_GasDisk_Dens[i][iri][izi] + tm1*App_GasDisk_Dens[i][iri-1][izi];
//         double dens2 = t*App_GasDisk_Dens[i][iri][izi-1] + tm1*App_GasDisk_Dens[i][iri-1][izi-1];
//         double dens3 = ss*App_GasDisk_Dens[i][iri][izi] + sm1*App_GasDisk_Dens[i][iri][izi-1];
//         double dens4 = ss*App_GasDisk_Dens[i][iri-1][izi] + sm1*App_GasDisk_Dens[i][iri-1][izi-1];
//     
//         //appdiskdens += 0.25*(dens1+dens2+dens3+dens4);
//         appdiskdens += ss*dens1+sm1*dens2;
//     }
//     
//     return appdiskdens;
}

//Get high frequency gas disk component potential (Psi^prime)
void GetAppGasDiskPot(void)
{
    cout << "Part 1" << endl;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        for (int j = 0; j < nr; ++j)
        {
            double temp_integral = 0;
            
            for (int k = 0; k < nr; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double dz = dr;
                double app_dens = GasDiskDensfI2(r, z, i);
                
                App_GasDisk_Dens[i][j][k] = app_dens;
                
                double integrand = app_dens*dz;
                temp_integral += integrand;
                
                App_GasDisk_Force[i][j][k] = temp_integral;
            }
        }
    }
    
    cout << "Part 1" << endl;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        for (int j = 0; j < nr; ++j)
        {
            double temp_integral = 0;
            
            for (int k = 0; k < nr; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double dz = dr;
                double app_force = App_GasDisk_Force[i][j][k];
                
                double integrand = app_force*dz;
                temp_integral += integrand;
                
                double rad = sqrt(r*r+z*z);
                double vtrunc = GetTruncGas(rad, i);
                
                App_GasDisk_Pot[i][j][k] = -temp_integral;//*vtrunc;
            }
        }
    }
    
    cout << "Part 1" << endl;
    
    //Now take gradient and Laplacian to get the high-frequency force and density
    for (int i = 0; i < gas_disk; ++i)
    {
        for (int j = 1; j < nr-1; ++j)
        {
            for (int k = 1; k < nr-1; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double delta = dr;
                double r1 = r+delta, r2 = r-delta;
                double z1 = z+delta, z2 = z-delta; 
                
                //cout << "here " << r << " " << z << endl;
                
                //Get dpsi/dr, dpsi/dz
                if (z != 0 && r != 0)
                {
                    double psir1 = App_GasDisk_Pot[i][j+1][k];
                    double psir2 = App_GasDisk_Pot[i][j-1][k];
                    double psiz1 = App_GasDisk_Pot[i][j][k+1];
                    double psiz2 = App_GasDisk_Pot[i][j][k-1];

                    DPsi_DR[i][j][k] = (psir1-psir2)/(r1-r2);
                    DPsi_DZ[i][j][k] = (psiz1-psiz2)/(z1-z2);

                    App_GasDisk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                    App_GasDisk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                }
                //else if (z == 0 && r != 0)
                //{
                //    DPsi_DR[i][j][k] = 0;
                //    DPsi_DZ[i][j][k] = 0;
                //    App_GasDisk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                //    App_GasDisk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                //}
                //else if (r == 0)//needs to be fixed, is not usually zero
                //{
                //    DPsi_DR[i][j][k] = 2*DPsi_DR[i][j][1]-DPsi_DR[i][j][2];
                //    DPsi_DZ[i][j][k] = 2*DPsi_DZ[i][j][1]-DPsi_DZ[i][j][2];
                //    App_GasDisk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                //    App_GasDisk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                //}
                //cout << "here " << j << " " << k << endl;
                
            }
        }
        
        //case z=0, z=r_max
        for (int j = 0; j < nr; ++j)
        {
            DPsi_DR[i][j][0] = 0;
            DPsi_DZ[i][j][0] = 0;
            App_GasDisk_Force_R[i][j][0] = DPsi_DR[i][j][0];
            App_GasDisk_Force_Z[i][j][0] = DPsi_DZ[i][j][0];
            
            DPsi_DR[i][j][nr-1] = 2*DPsi_DR[i][j][nr-2]-DPsi_DR[i][j][nr-3];
            DPsi_DZ[i][j][nr-1] = 2*DPsi_DZ[i][j][nr-2]-DPsi_DZ[i][j][nr-3];
            App_GasDisk_Force_R[i][j][nr-1] = DPsi_DR[i][j][nr-1];
            App_GasDisk_Force_Z[i][j][nr-1] = DPsi_DZ[i][j][nr-1];
        }
        
        //case R=0, R=r_max
        for (int k = 0; k < nr; ++k)
        {
            DPsi_DR[i][0][k] = 2*DPsi_DR[i][1][k]-DPsi_DR[i][2][k];
            DPsi_DZ[i][0][k] = 2*DPsi_DZ[i][1][k]-DPsi_DZ[i][2][k];
            App_GasDisk_Force_R[i][0][k] = DPsi_DR[i][0][k];
            App_GasDisk_Force_Z[i][0][k] = DPsi_DZ[i][0][k];
            
            DPsi_DR[i][nr-1][k] = 2*DPsi_DR[i][nr-2][k]-DPsi_DR[i][nr-3][k];
            DPsi_DZ[i][nr-1][k] = 2*DPsi_DZ[i][nr-2][k]-DPsi_DZ[i][nr-3][k];
            App_GasDisk_Force_R[i][nr-1][k] = DPsi_DR[i][nr-1][k];
            App_GasDisk_Force_Z[i][nr-1][k] = DPsi_DZ[i][nr-1][k];
        }
    }
    
    cout << "Part 2" << endl;
    
    for (int i = 0; i < gas_disk; ++i)
    {
        for (int j = 1; j < nr-1; ++j)
        {
            for (int k = 1; k < nr-1; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double delta = dr;
                double r1 = r+delta, r2 = r-delta;
                double z1 = z+delta, z2 = z-delta; 
                
                //Get d2psi/dr2, d2psi/dz2
                if (z != 0 && r != 0)
                {
                    double dpsidr1 = DPsi_DR[i][j+1][k];
                    double dpsidr2 = DPsi_DR[i][j-1][k];
                    double dpsidz1 = DPsi_DZ[i][j][k+1];
                    double dpsidz2 = DPsi_DZ[i][j][k-1];

                    D2Psi_DR2[i][j][k] = (dpsidr1-dpsidr2)/(r1-r2);
                    D2Psi_DZ2[i][j][k] = (dpsidz1-dpsidz2)/(z1-z2);

                    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                                       D2Psi_DZ2[i][j][k];

                    App_GasDisk_Dens[i][j][k] = -laplacian*oneover4pi;
                }
                //else if (z == 0)
                //{
                //    D2Psi_DR2[i][j][k] = 2*D2Psi_DR2[i][1][k]-D2Psi_DR2[i][2][k];
                //    D2Psi_DZ2[i][j][k] = 2*D2Psi_DZ2[i][1][k]-D2Psi_DZ2[i][2][k];
                //    
                //    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                //                       D2Psi_DZ2[i][j][k];
                // 
                //    App_GasDisk_Dens[i][j][k] = -laplacian*oneover4pi;
                //}
                //else if (r == 0)//needs to be fixed, is not usually zero
                //{
                //    D2Psi_DR2[i][j][k] = 2*D2Psi_DR2[i][j][1]-D2Psi_DR2[i][j][2];
                //    D2Psi_DZ2[i][j][k] = 2*D2Psi_DZ2[i][j][1]-D2Psi_DZ2[i][j][2];
                //    
                //    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                //                       D2Psi_DZ2[i][j][k];
                //
                //    App_GasDisk_Dens[i][j][k] = -laplacian*oneover4pi;
                //}
            }
        }
        
        //case z=0, z=r_max
        for (int j = 0; j < nr; ++j)
        {
            D2Psi_DR2[i][j][0] = 2*D2Psi_DR2[i][j][1]-D2Psi_DR2[i][j][2];
            D2Psi_DZ2[i][j][0] = 2*D2Psi_DZ2[i][j][1]-D2Psi_DZ2[i][j][2];
            
            double laplacian = D2Psi_DR2[i][j][0] + //DPsi_DR[i][j][0]/r + 
                               D2Psi_DZ2[i][j][0];

            App_GasDisk_Dens[i][j][0] = -laplacian*oneover4pi;
            
            D2Psi_DR2[i][j][nr-1] = 2*D2Psi_DR2[i][j][nr-2]-D2Psi_DR2[i][j][nr-3];
            D2Psi_DZ2[i][j][nr-1] = 2*D2Psi_DZ2[i][j][nr-2]-D2Psi_DZ2[i][j][nr-3];
            
            laplacian = D2Psi_DR2[i][j][nr-1] + DPsi_DR[i][j][nr-1]/Radius[nr-1] + 
                        D2Psi_DZ2[i][j][nr-1];

            App_GasDisk_Dens[i][j][nr-1] = -laplacian*oneover4pi;
        }
        
        //case R=0, R=r_max
        for (int k = 0; k < nr; ++k)
        {
            D2Psi_DR2[i][0][k] = 2*D2Psi_DR2[i][1][k]-D2Psi_DR2[i][2][k];
            D2Psi_DZ2[i][0][k] = 2*D2Psi_DZ2[i][1][k]-D2Psi_DZ2[i][2][k];
            
            double laplacian = D2Psi_DR2[i][0][k] + //DPsi_DR[i][0][k]/r + 
                               D2Psi_DZ2[i][0][k];

            App_GasDisk_Dens[i][0][k] = -laplacian*oneover4pi;
            
            D2Psi_DR2[i][nr-1][k] = 2*D2Psi_DR2[i][nr-2][k]-D2Psi_DR2[i][nr-3][k];
            D2Psi_DZ2[i][nr-1][k] = 2*D2Psi_DZ2[i][nr-2][k]-D2Psi_DZ2[i][nr-3][k];
            
            laplacian = D2Psi_DR2[i][nr-1][k] + DPsi_DR[i][nr-1][k]/Radius[nr-1] + 
                        D2Psi_DZ2[i][nr-1][k];

            App_GasDisk_Dens[i][nr-1][k] = -laplacian*oneover4pi;
        }
    }
}
    
                
    
    
void GetAppDiskPot(void)
{
    cout << "Part 1" << endl;
    
    for (int i = 0; i < disk; ++i)
    {
        for (int j = 0; j < nr; ++j)
        {
            double temp_integral = 0;
            
            for (int k = 0; k < nr; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double dz = dr;
                double app_dens = DiskDensfI(r, z, i);
                
                App_Disk_Dens[i][j][k] = app_dens;
                
                double integrand = app_dens*dz;
                temp_integral += integrand;
                
                App_Disk_Force[i][j][k] = temp_integral;
                if (k==0)
                {
                    App_Disk_Force[i][j][k] = 0;
                }
            }
        }
    }
    
    cout << "Part 1" << endl;
    
    for (int i = 0; i < disk; ++i)
    {
        for (int j = 0; j < nr; ++j)
        {
            double temp_integral = 0;
            
            for (int k = 0; k < nr; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double dz = dr;
                double app_force = App_Disk_Force[i][j][k];
                
                double integrand = app_force*dz;
                double rad = sqrt(r*r+z*z);
                double vtrunc = GetTrunc(rad, i);
                
                temp_integral += integrand;
                
                App_Disk_Pot[i][j][k] = -fourpi*temp_integral;//*vtrunc;
                if (k==0)
                {
                    App_Disk_Pot[i][j][k] = 0;
                }
                //if(k<30)
                //    cout << i << " " << r << " " << z << " " 
                //         << App_Disk_Dens[i][j][k] << " " << AppDiskDens3(r,z) << "     "
                //         << App_Disk_Pot[i][j][k] << " " << AppDiskPot3(r,z) << " " 
                //         << 1.68*log(cosh(z/G.Z_Disk[i]))*exp(r/G.R_Disk[i]) << endl;
            }
        }
    }
    
    cout << "Part 1" << endl;
    
    //Now take gradient and Laplacian to get the high-frequency force and density
    for (int i = 0; i < disk; ++i)
    {
        for (int j = 1; j < nr-1; ++j)
        {
            for (int k = 1; k < nr-1; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double delta = dr;
                double r1 = r+delta, r2 = r-delta;
                double z1 = z+delta, z2 = z-delta; 
                
                //cout << "here " << r << " " << z << endl;
                
                //Get dpsi/dr, dpsi/dz
                if (z != 0 && r != 0)
                {
                    double psir1 = App_Disk_Pot[i][j+1][k];
                    double psir2 = App_Disk_Pot[i][j-1][k];
                    double psiz1 = App_Disk_Pot[i][j][k+1];
                    double psiz2 = App_Disk_Pot[i][j][k-1];

                    DPsi_DR[i][j][k] = (psir1-psir2)/(r1-r2);
                    DPsi_DZ[i][j][k] = (psiz1-psiz2)/(z1-z2);

                    App_Disk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                    App_Disk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                }
                //else if (z == 0 && r != 0)
                //{
                //    DPsi_DR[i][j][k] = 0;
                //    DPsi_DZ[i][j][k] = 0;
                //    App_Disk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                //    App_Disk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                //}
                //else if (r == 0)//needs to be fixed, is not usually zero
                //{
                //    DPsi_DR[i][j][k] = 2*DPsi_DR[i][j][1]-DPsi_DR[i][j][2];
                //    DPsi_DZ[i][j][k] = 2*DPsi_DZ[i][j][1]-DPsi_DZ[i][j][2];
                //    App_Disk_Force_R[i][j][k] = DPsi_DR[i][j][k];
                //    App_Disk_Force_Z[i][j][k] = DPsi_DZ[i][j][k];
                //}
                //cout << "here " << j << " " << k << endl;
                
            }
        }
        
        //case z=0, z=r_max
        for (int j = 0; j < nr; ++j)
        {
            DPsi_DR[i][j][0] = 0;
            DPsi_DZ[i][j][0] = 0;
            App_Disk_Force_R[i][j][0] = DPsi_DR[i][j][0];
            App_Disk_Force_Z[i][j][0] = DPsi_DZ[i][j][0];
            
            DPsi_DR[i][j][nr-1] = 2*DPsi_DR[i][j][nr-2]-DPsi_DR[i][j][nr-3];
            DPsi_DZ[i][j][nr-1] = 2*DPsi_DZ[i][j][nr-2]-DPsi_DZ[i][j][nr-3];
            App_Disk_Force_R[i][j][nr-1] = DPsi_DR[i][j][nr-1];
            App_Disk_Force_Z[i][j][nr-1] = DPsi_DZ[i][j][nr-1];
        }
        
        //case R=0, R=r_max
        for (int k = 0; k < nr; ++k)
        {
            DPsi_DR[i][0][k] = 2*DPsi_DR[i][1][k]-DPsi_DR[i][2][k];
            DPsi_DZ[i][0][k] = 2*DPsi_DZ[i][1][k]-DPsi_DZ[i][2][k];
            App_Disk_Force_R[i][0][k] = DPsi_DR[i][0][k];
            App_Disk_Force_Z[i][0][k] = DPsi_DZ[i][0][k];
            
            DPsi_DR[i][nr-1][k] = 2*DPsi_DR[i][nr-2][k]-DPsi_DR[i][nr-3][k];
            DPsi_DZ[i][nr-1][k] = 2*DPsi_DZ[i][nr-2][k]-DPsi_DZ[i][nr-3][k];
            App_Disk_Force_R[i][nr-1][k] = DPsi_DR[i][nr-1][k];
            App_Disk_Force_Z[i][nr-1][k] = DPsi_DZ[i][nr-1][k];
        }
    }
    
    cout << "Part 2" << endl;
    
    for (int i = 0; i < disk; ++i)
    {
        for (int j = 1; j < nr-1; ++j)
        {
            for (int k = 1; k < nr-1; ++k)
            {
                double r = Radius[j];
                double z = Radius[k];
                double delta = dr;
                double r1 = r+delta, r2 = r-delta;
                double z1 = z+delta, z2 = z-delta; 
                
                //Get d2psi/dr2, d2psi/dz2
                if (z != 0 && r != 0)
                {
                    double dpsidr1 = DPsi_DR[i][j+1][k];
                    double dpsidr2 = DPsi_DR[i][j-1][k];
                    double dpsidz1 = DPsi_DZ[i][j][k+1];
                    double dpsidz2 = DPsi_DZ[i][j][k-1];

                    D2Psi_DR2[i][j][k] = (dpsidr1-dpsidr2)/(r1-r2);
                    D2Psi_DZ2[i][j][k] = (dpsidz1-dpsidz2)/(z1-z2);

                    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                                       D2Psi_DZ2[i][j][k];

                    App_Disk_Dens[i][j][k] = -laplacian*oneover4pi;
                    
                    //if(k<20)
                    //    cout << i << " " << r << " " << z << " " 
                    //         << App_Disk_Dens[i][j][k] << " " << AppDiskDens3(r,z) << "     "
                    //         << App_Disk_Pot[i][j][k] << " " << AppDiskPot3(r,z) << endl;
                }
                //else if (z == 0)
                //{
                //    D2Psi_DR2[i][j][k] = 2*D2Psi_DR2[i][1][k]-D2Psi_DR2[i][2][k];
                //    D2Psi_DZ2[i][j][k] = 2*D2Psi_DZ2[i][1][k]-D2Psi_DZ2[i][2][k];
                //    
                //    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                //                       D2Psi_DZ2[i][j][k];
                // 
                //    App_Disk_Dens[i][j][k] = -laplacian*oneover4pi;
                //}
                //else if (r == 0)//needs to be fixed, is not usually zero
                //{
                //    D2Psi_DR2[i][j][k] = 2*D2Psi_DR2[i][j][1]-D2Psi_DR2[i][j][2];
                //    D2Psi_DZ2[i][j][k] = 2*D2Psi_DZ2[i][j][1]-D2Psi_DZ2[i][j][2];
                //    
                //    double laplacian = D2Psi_DR2[i][j][k] + DPsi_DR[i][j][k]/r + 
                //                       D2Psi_DZ2[i][j][k];
                //
                //    App_Disk_Dens[i][j][k] = -laplacian*oneover4pi;
                //}
            }
        }
        
        //case z=0, z=r_max
        for (int j = 0; j < nr; ++j)
        {
            D2Psi_DR2[i][j][0] = 2*D2Psi_DR2[i][j][1]-D2Psi_DR2[i][j][2];
            D2Psi_DZ2[i][j][0] = 2*D2Psi_DZ2[i][j][1]-D2Psi_DZ2[i][j][2];
            
            double laplacian = D2Psi_DR2[i][j][0] + //DPsi_DR[i][j][0]/r + 
                               D2Psi_DZ2[i][j][0];

            App_Disk_Dens[i][j][0] = laplacian*oneover4pi;
            
            D2Psi_DR2[i][j][nr-1] = 2*D2Psi_DR2[i][j][nr-2]-D2Psi_DR2[i][j][nr-3];
            D2Psi_DZ2[i][j][nr-1] = 2*D2Psi_DZ2[i][j][nr-2]-D2Psi_DZ2[i][j][nr-3];
            
            laplacian = D2Psi_DR2[i][j][nr-1] + DPsi_DR[i][j][nr-1]/Radius[nr-1] + 
                        D2Psi_DZ2[i][j][nr-1];

            App_Disk_Dens[i][j][nr-1] = laplacian*oneover4pi;
        }
        
        //case R=0, R=r_max
        for (int k = 0; k < nr; ++k)
        {
            D2Psi_DR2[i][0][k] = 2*D2Psi_DR2[i][1][k]-D2Psi_DR2[i][2][k];
            D2Psi_DZ2[i][0][k] = 2*D2Psi_DZ2[i][1][k]-D2Psi_DZ2[i][2][k];
            
            double laplacian = D2Psi_DR2[i][0][k] + //DPsi_DR[i][0][k]/r + 
                               D2Psi_DZ2[i][0][k];

            App_Disk_Dens[i][0][k] = laplacian*oneover4pi;
            
            D2Psi_DR2[i][nr-1][k] = 2*D2Psi_DR2[i][nr-2][k]-D2Psi_DR2[i][nr-3][k];
            D2Psi_DZ2[i][nr-1][k] = 2*D2Psi_DZ2[i][nr-2][k]-D2Psi_DZ2[i][nr-3][k];
            
            laplacian = D2Psi_DR2[i][nr-1][k] + DPsi_DR[i][nr-1][k]/Radius[nr-1] + 
                        D2Psi_DZ2[i][nr-1][k];

            App_Disk_Dens[i][nr-1][k] = laplacian*oneover4pi;
        }
    }
}
    
                
    
    
