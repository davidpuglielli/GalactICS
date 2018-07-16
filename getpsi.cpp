//The iniital versions of these functions included redundant calculations, eg.
//getting t three times for three components with the same radius. Now we've
//merged them into one function to avoid that extra overhead.

#include "galaxy.h"

double GetTotalPsi(double &r)
{
    double gettotalpsi = 0;
    
    double log_r = log10(r);
    //int ihi = ceil((log_r-log_dr)/delta_logr+1);
    int ihi = ceil(r/dr);
    //cout <<setprecision(8)<<"ihi "<<ihi<<" "<<r<<" "<<int((log_r-log_dr)/delta_logr)+1<<" "<<log_rmax << endl;
    
    //if (ihi < 0) ihi = 0;
    //else if (ihi > nr-1) ihi = nr-1;
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetTotalPsi finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
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
    double tm = 1 - t;
    
    //cout <<"ihi "<<ihi<<" "<<r<< " "<<r1<<" "<<r2<<" "<<t<<" "<<tm<<endl;
    
    if (disk_flag == 1)
    {
        gettotalpsi += t*D_Pot[ihi] + tm*D_Pot[ihi-1];
    }
    if (gasdisk_flag == 1)
    {
        gettotalpsi += t*G_Pot[ihi] + tm*G_Pot[ihi-1];
    }
    if (bulge_flag == 1)
    {
        if (sersic_flag)
        {
            gettotalpsi += SersicPotential(r);
        }
        else
        {
            gettotalpsi += t*B_Pot[ihi] + tm*B_Pot[ihi-1];
        }
    }
    if (halo_flag == 1)
    {
        gettotalpsi += t*H_Pot[ihi] + tm*H_Pot[ihi-1];
    }
    
    return gettotalpsi;
}

double GetDiskPsi(double &r, double &z)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetDiskPsi finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
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
    double tm = 1 - t;
	
    return t*D_Pot[ihi] + tm*D_Pot[ihi-1];
}

double GetGasDiskPsi(double &r, double &z)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetDiskPsi finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
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
    double tm = 1 - t;
	
    return t*G_Pot[ihi] + tm*G_Pot[ihi-1];
}
double GetHaloPsi(double &r)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetHaloPsi finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
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
    double tm = 1 - t;
	
    return t*H_Pot[ihi] + tm*H_Pot[ihi-1];
}


double GetBulgePsi(double &r)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GetBulgePsi finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
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
    double tm = 1 - t;
	
    return t*B_Pot[ihi] + tm*B_Pot[ihi-1];
}

        
