#include "galaxy.h"

void HaloPotentialEstimate(void)
{
    double eps = 0.0001;
    
    for (int i = 0; i < nr; ++i)
    {
        double r = Radius[i];
        if (r == 0)
            r = eps;
        H_Dens[i]=HaloDensity(r);
    }
    
// KK: now get the potential harmonics of this new density. (BT 2-208)
// Simpson's rule integration.

    double ds = delta_logr;
    
    S_1[0] = 0;
    double r = Radius[2];
    double ss = log10(r);
    
    S_1[2] = (r*dr/3)*(4*H_Dens[1]*(1-dr/r)*(1-dr/r)+H_Dens[2]);
    //S_1[2] = (ds/3) * (4*H_Dens[1]*exp(3*ss-3*ds) + H_Dens[2]*exp(3*ss));cout<<"s1   "<<S_1[2]<<endl;
    
    double r_old = r;
    
    for (int ir = 4; ir < nr; ir+=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s1a = (ds/3)*(H_Dens[ir-2]*exp(3*s-6*ds) + 4*H_Dens[ir-1]*exp(3*s-3*ds) +
        //                     H_Dens[ir]*exp(3*s));
        double s1a = (r*dr/3)*(H_Dens[ir-2]*pow(1-2*dr/r, 2) + 
                     4*H_Dens[ir-1]*pow(1-dr/r, 2) + H_Dens[ir]);
        S_1[ir] = s1a+S_1[ir-2]*r_old/r;
        r_old = r;//cout<<"s1   "<<r<<" "<<s1a<<" "<<S_1[ir]<<endl;
    }
    
    S_2[nr-1] = 0;
    
    r_old = Radius[nr-1];
    
    for (int ir = nr-3; ir > 1; ir-=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s2a = (ds/3)*(H_Dens[ir]*exp(3*s) + 4*H_Dens[ir+1]*exp(3*s+3*ds) +
        //                     H_Dens[ir+2]*exp(3*s+6*ds));
        double s2a = (r*dr/3) * (H_Dens[ir+2]*(1+2*dr/r) + 
                     4*H_Dens[ir+1]*(1+dr/r) + H_Dens[ir]);
        S_2[ir] = s2a + S_2[ir+2];//cout<<"s2 h   "<<r<<" "<<s2a<<" "<<S_2[ir]<<endl;
        r_old = r;
    }
    
    for (int ir = 2; ir < nr; ir+=2)
    {
        r = Radius[ir];
        H_Pot[ir] = fourpi*(S_1[ir]+S_2[ir]);
        H_FR[ir] = -fourpi*S_1[ir]/r;
#ifdef DEBUG
        //cout << "Halopotentialestimate  " << r << " " << ir << " " << H_Pot[ir] << " " << H_FR[ir] << endl;
#endif
    }
    
    H_Pot[0] = 3*(H_Pot[2]-H_Pot[4])+H_Pot[6];
    //cout << "Halopotentialestimate  " << r << " " << H_Pot[0] << " " << H_FR[0] << endl;
    H_FR[0] = 0;
    
//  KK: then linearly interpolate other bins.
    
    for (int ir = 1; ir < nr; ir+=2)
    {
        H_Pot[ir] = (H_Pot[ir-1]+H_Pot[ir+1])/2;
        H_FR[ir] = (H_FR[ir-1]+H_FR[ir+1])/2;
#ifdef DEBUG
        //cout << "Halopotentialestimate  " << Radius[ir] << " " << ir << " " << H_Pot[ir] << " " << H_FR[ir] << endl;
#endif
    }
}

