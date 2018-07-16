#include "galaxy.h"
            
double HaloForce(double &r)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/delta_logr+1);// /(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "GettotalPsi finds out of range indices. Exiting..." << endl;
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
    
    return t*H_FR[ihi] + tm1*H_FR[ihi-1];
}

            
            
            
            
    
        
