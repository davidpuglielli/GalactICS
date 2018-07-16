//KK:Returns the potential, being the sum of the spherical harmonics
//with coefficients in apot, and the potential appdiskpot which
//approximates the high-frequency components of the disk.

#include "galaxy.h"
#include <limits>

double Pot(double s, double z)
{
    double pot = 0, r=sqrt(s*s+z*z);
    
    if (r == 0)
    {
        return pot = A_Pot[0][0]*oneoversqrt4pi;
    }
    
    int ihi = ceil(r/dr);
	
    if(r < dr)
    {
        ihi = 1;
    }
    else if (r > 10*r_max)
    {
        //cout << "radius too high in Pot. Continuing... " << r << " " 
        //     << s << " "<< z << endl;
        ihi = nr-1;
    }
    else if (ihi < 1)
    {
        cout << "Pot finds out of range indices. Possible int overflow. " <<
                "Exiting..." << endl;
        cout <<"ihi = "<<ihi<<" "<<r<<" "<<dr<<" "<<r/dr<<" "<<ceil(r/dr)<<" "
             <<s<<" "<<z<<endl;
        exit(1);
    }
    else if (ihi > nr-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<s<<" "<<z<<" "<< endl;
        ihi = nr-1;
    }
    
    double r1 = Radius[ihi-1];
    double r2 = Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    double cos_theta = z/r;
    
    for (int l = l_max; l > -1; l-=2)
    {
        //pot += Legendre_Table[l/2][int(1000*cos_theta)]*Plcon[l]*
        //       (t*A_Pot[l/2][ihi] + tm1*A_Pot[l/2][ihi-1]);
        pot += gsl_sf_legendre_Pl(l, cos_theta)*Plcon[l]*
               (t*A_Pot[l/2][ihi] + tm1*A_Pot[l/2][ihi-1]);
        
        //if(nbody_flag==1)       
        //if(z<1&&r<1&&iter==0)
        //cout<<"pot   "<<s<<" "<<z<<" "<<l<<" "<<cos_theta<<" "
        //    <<gsl_sf_legendre_Pl(l, cos_theta)<<" "<<Plcon[l]<<"\n"
        //    <<"      "<<ihi<<" "<<t<<" "<<tm1<<" "<<A_Pot[l/2][ihi]<<" "<<A_Pot[l/2][ihi-1]<<endl;
    }
    
    if (disk_flag)
    {
        pot += AppDiskPot(s, z);
    }
    
    if (gasdisk_flag)
    {
	    pot += AppGasDiskPot(s, z);
    }
    
    if (smbh_flag)
    {
        pot += G.bh_mass/sqrt(s*s+z*z);
    }
    
    return pot;
}

//Do the same thing as above but for any array of harmonics. Used in the GetFreqs function
double PotVector(double s, double z, vector<vector<double> >& Vec)
{
    //double linearp, quadraticp, rfit[10], pfit[10];
    
    double pot = 0, r=sqrt(s*s+z*z);
    
    if (r == 0)
    {
        return pot = Vec[0][0]*oneoversqrt4pi;
    }
    
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Pot finds out of range indices. Exiting..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        exit(1);
    }
    else if (ihi > nr-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr-1;
    }
    
    double r1 = Radius[ihi-1];
    double r2 = Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    double cos_theta = z/r;
    
    //cout << "PotVec   ";
    
    for (int l = l_max; l > -1; l-=2)
    {
        pot += gsl_sf_legendre_Pl(l, cos_theta)*Plcon[l]*
               (t*Vec[l/2][ihi] + tm1*Vec[l/2][ihi-1]);
    }
    
    if (disk_flag == 1)
    {
        pot += AppDiskPot(s, z);
    }
    
    if (gasdisk_flag)
    {
	    pot += AppGasDiskPot(s, z);
    }
    
    return pot;
}

