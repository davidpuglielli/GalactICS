//NB. Need to use nr+1 in for loops because it solves issues regarding the 
//interpolation with an even number of bins

#include "galaxy.h"

void DiskPotentialEstimate(void)
{
    double s = 0, eps = 0.1*dr;
    int n_theta = 100;
    
    //DiskPotentialEstimateCUDA();
    
    for (int i = 0; i < nr; ++i)
    {
        double r = Radius[i];
        //double s=0;
        
        if (r == 0)
            r = eps;
        
        double cos_theta_max = min(1.0, 10*G.Z_Disk[0]/r);
        double d_cos_theta = cos_theta_max / n_theta;
        
        s+= dPolarDiskDens(r, cos_theta_max) + dPolarDiskDens(r, 0);
        
        for (int is = 1; is < n_theta; is+=2)
        {
            double cos_theta = is*d_cos_theta;
            s += 4*dPolarDiskDens(r, cos_theta);
        }
        for (int is = 2; is < n_theta-1; is+=2)
        {
            double cos_theta = is*d_cos_theta;
            s += 2*dPolarDiskDens(r, cos_theta);
        }
        
        s *= d_cos_theta/3;
        D_Dens[i] = s;
    }
    
    //KK: now get the potential harmonics of this new density. (BT 2-208)
    //Simpson's rule integration.

    double ds = delta_logr;
    
    S_1[0] = 0;
    double r = Radius[2];
    double ss = log10(r);
    S_1[2] = (r*dr/3)*(4*D_Dens[1]*(1-dr/r)*(1-dr/r)+D_Dens[2]);
    //S_1[2] = (ds/3) * (4*D_Dens[1]*exp(3*ss-3*ds) + D_Dens[2]*exp(3*ss));cout<<"s1   "<<S_1[2]<<endl;
    
    double r_old = r;
    
    for (int ir = 4; ir < nr; ir+=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s1a = (ds/3)*(D_Dens[ir-2]*exp(3*s-6*ds) + 4*D_Dens[ir-1]*exp(3*s-3*ds) +
        //                     D_Dens[ir]*exp(3*s));
        double s1a = (r*dr/3)*(D_Dens[ir-2]*pow(1-2*dr/r, 2) + 
                     4*D_Dens[ir-1]*pow(1-dr/r, 2) + D_Dens[ir]);
        S_1[ir] = s1a+S_1[ir-2]*r_old/r;//cout<<"s1 d   "<<r<<" "<<s1a<<" "<<S_1[ir]<<" "<<D_Dens[ir-2]<<endl;
        r_old = r;
    }
    
    S_2[nr-1] = 0;
    
    r_old = Radius[nr-1];
    
    for (int ir = nr-3; ir > 1; ir-=2)//NB. We are hitting even numbers in this loop
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s2a = (ds/3)*(D_Dens[ir]*exp(3*s) + 4*D_Dens[ir+1]*exp(3*s+3*ds) +
        //                     D_Dens[ir+2]*exp(3*s+6*ds));
        double s2a = (r*dr/3)*(D_Dens[ir+2]*(1+2*dr/r) + 
                     4*D_Dens[ir+1]*(1+dr/r) + D_Dens[ir]);
        S_2[ir] = s2a+S_2[ir+2];//cout<<"s2 d   "<<r<<" "<<s2a<<" "<<S_2[ir]<<" "<<D_Dens[ir+2]<<endl;
        r_old = r;
    }
    
    for (int ir = 2; ir < nr; ir+=2)
    {
        r = Radius[ir];
        D_Pot[ir] = fourpi*(S_1[ir]+S_2[ir]);
        D_FR[ir] = -fourpi*S_1[ir]/r;
#ifdef DEBUG
        //cout << "Diskpotentialestimate  " << r << " " << ir << " " << D_Pot[ir] << " " << D_FR[ir] << endl;
#endif
    }
    
    D_Pot[0] = 3*(D_Pot[2]-D_Pot[4])+D_Pot[6];
    D_FR[0] = 0;
    //cout<<D_Pot[0]<<" "<<D_FR[0]<<endl;
    
//  KK: then linearly interpolate other bins.
    
    for (int ir = 1; ir < nr; ir+=2)
    {
        D_Pot[ir] = (D_Pot[ir-1]+D_Pot[ir+1])/2;
        D_FR[ir] = (D_FR[ir-1]+D_FR[ir+1])/2;
#ifdef DEBUG
        //cout << "Diskpotest " << Radius[ir] << " " << ir << " " << D_Pot[ir] << " " << D_FR[ir] << endl;
#endif
    }
}
    
double DiskForce(double &r)
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
    
    //double cos_theta = z/r;
    //l_max = lmax;

    return t*D_FR[ihi] + tm1*D_FR[ihi-1];
}

double DiskDensity(double &r)
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
    
    //double cos_theta = z/r;
    //l_max = lmax;

    return t*D_Dens[ihi] + tm1*D_Dens[ihi-1];
}
   
double GetDiskMass(int i)
{
    //Get the mass of each disk using Simpson's rule
    if (n_simpson%2==1)
        ++n_simpson;
    
    vector<double> Annulus_Mass(n_simpson, 0);
    double disk_mass = 0;
    double drad = Disk_Edge[i]/n_simpson, dz = 10*G.Z_Disk[i]/n_simpson;
    
    for (int j = 0; j < n_simpson; ++j)
    {
        double radius = j*drad;
        
        Annulus_Mass.at(j) += DiskDensfI(radius, 0, i)+4*DiskDensfI(radius, dz, i)+
                              DiskDensfI(radius, 5*G.Z_Disk[i], i);

        for (int k = 2; k < n_simpson; k+=2)
        {
            double z = k*dz;
            Annulus_Mass.at(j) += 2*DiskDensfI(radius, z, i)+4*DiskDensfI(radius, z+dz, i);
        }
        
        Annulus_Mass.at(j) *= dz/3;
        Annulus_Mass.at(j) *= radius;
    }
        
    //Now apply Simpson's rule in the radial direction
    disk_mass += Annulus_Mass.at(0)+4*Annulus_Mass.at(1)+Annulus_Mass.at(n_simpson-1);

    for (int j = 2; j < n_simpson; j+=2)
    {
        disk_mass += 2*Annulus_Mass.at(j)+4*Annulus_Mass.at(j+1);
    }
    
    disk_mass *= drad*fourpi/3;
    
    return disk_mass;
}


void GasDiskPotentialEstimate(void)
{
    double s = 0, eps = 0.1*dr;
    int n_theta = 100;
    
    //cout << "Gas " << endl;
    
    for (int i = 0; i < nr; ++i)
    {
        //cout << "Gas " << nr << endl;
        double r = Radius[i];
        
        if (r == 0)
            r = eps;
        
        double cos_theta_max = min(1.0, 10*G.Z_GasDisk[0]/r);
        double d_cos_theta = cos_theta_max / n_theta;
        
        s+= dPolarGasDiskDens(r, cos_theta_max) + dPolarGasDiskDens(r, 0);
        
        //cout << "Gas " << dPolarGasDiskDens(r, cos_theta_max) << endl;
        for (int is = 1; is < n_theta; is+=2)
        {
            double cos_theta = is*d_cos_theta;
            s += 4*dPolarGasDiskDens(r, cos_theta);
            //cout << "Gas " << is << " " << dPolarGasDiskDens(r, cos_theta_max) << endl;
        }
        for (int is = 2; is < n_theta-1; is+=2)
        {
            double cos_theta = is*d_cos_theta;
            s += 2*dPolarGasDiskDens(r, cos_theta);
            //cout << "Gas " << is << " " << dPolarGasDiskDens(r, cos_theta_max) << endl;
        }
        
        s *= d_cos_theta/3;
        G_Dens[i] = s;
    }
    
    //KK: now get the potential harmonics of this new density. (BT 2-208)
    //Simpson's rule integration.

    double ds = delta_logr;
    
    S_1[0] = 0;
    double r = Radius[2];
    double ss = log10(r);
    S_1[2] = (r*dr/3)*(4*G_Dens[1]*(1-dr/r)*(1-dr/r)+G_Dens[2]);
    //S_1[2] = (ds/3) * (4*D_Dens[1]*exp(3*ss-3*ds) + D_Dens[2]*exp(3*ss));cout<<"s1   "<<S_1[2]<<endl;
    
    double r_old = r;
    
    for (int ir = 4; ir < nr; ir+=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s1a = (ds/3)*(D_Dens[ir-2]*exp(3*s-6*ds) + 4*D_Dens[ir-1]*exp(3*s-3*ds) +
        //                     D_Dens[ir]*exp(3*s));
        double s1a = (r*dr/3)*(G_Dens[ir-2]*pow(1-2*dr/r, 2) + 
                     4*G_Dens[ir-1]*pow(1-dr/r, 2) + G_Dens[ir]);
        S_1[ir] = s1a+S_1[ir-2]*r_old/r;//cout<<"s1 d   "<<r<<" "<<s1a<<" "<<S_1[ir]<<" "<<D_Dens[ir-2]<<endl;
        r_old = r;
    }
    
    S_2[nr-1] = 0;
    
    r_old = Radius[nr-1];
    
    for (int ir = nr-3; ir > 1; ir-=2)//NB. We are hitting even numbers in this loop
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s2a = (ds/3)*(D_Dens[ir]*exp(3*s) + 4*D_Dens[ir+1]*exp(3*s+3*ds) +
        //                     D_Dens[ir+2]*exp(3*s+6*ds));
        double s2a = (r*dr/3)*(G_Dens[ir+2]*(1+2*dr/r) + 
                     4*G_Dens[ir+1]*(1+dr/r) + G_Dens[ir]);
        S_2[ir] = s2a+S_2[ir+2];//cout<<"s2 d   "<<r<<" "<<s2a<<" "<<S_2[ir]<<" "<<D_Dens[ir+2]<<endl;
        r_old = r;
    }
    
    for (int ir = 2; ir < nr; ir+=2)
    {
        r = Radius[ir];
        G_Pot[ir] = fourpi*(S_1[ir]+S_2[ir]);
        G_FR[ir] = -fourpi*S_1[ir]/r;
#ifdef DEBUG
        //cout << "Diskpotentialestimate  " << r << " " << ir << " " << D_Pot[ir] << " " << D_FR[ir] << endl;
#endif
    }
    
    G_Pot[0] = 3*(G_Pot[2]-G_Pot[4])+G_Pot[6];
    G_FR[0] = 0;
    //cout<<D_Pot[0]<<" "<<D_FR[0]<<endl;
    
//  KK: then linearly interpolate other bins.
    
    for (int ir = 1; ir < nr; ir+=2)
    {
        G_Pot[ir] = (G_Pot[ir-1]+G_Pot[ir+1])/2;
        G_FR[ir] = (G_FR[ir-1]+G_FR[ir+1])/2;
    }
}
    
double GasDiskForce(double &r)
{
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    
    if (ihi < 0) ihi = 0;
    else if (ihi > nr-1) ihi = nr-1;
    
    double r1 = Radius[ihi-1];
    double r2 = Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    //double cos_theta = z/r;
    //l_max = lmax;

    return t*G_FR[ihi] + tm1*G_FR[ihi-1];
}

double GasDiskDensity(double &r)
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
    
    //double cos_theta = z/r;
    //l_max = lmax;

    return t*G_Dens[ihi] + tm1*G_Dens[ihi-1];
}
   
double GetGasDiskMass(int i)
{
    //Get the mass of each disk using Simpson's rule
    if (n_simpson%2==1)
        ++n_simpson;
    
    vector<double> Annulus_Mass(n_simpson, 0);
    double disk_mass = 0;
    double drad = GasDisk_Edge[i]/n_simpson, dz = 50*G.Z_GasDisk[i]/n_simpson;
    
    for (int j = 0; j < n_simpson; ++j)
    {
        double radius = j*drad;
        
        Annulus_Mass.at(j) += GasDiskDensfI2(radius, 0, i)+
                              4*GasDiskDensfI2(radius, dz, i)+
                              GasDiskDensfI2(radius, 25*G.Z_GasDisk[i], i);

        for (int k = 2; k < n_simpson; k+=2)
        {
            double z = k*dz;
            Annulus_Mass.at(j) += 2*GasDiskDensfI2(radius, z, i)+
                                  4*GasDiskDensfI2(radius, z+dz, i);
        }
        
        Annulus_Mass.at(j) *= dz/3;
        Annulus_Mass.at(j) *= radius;
    }
        
    //Now apply Simpson's rule in the radial direction
    disk_mass += Annulus_Mass.at(0)+4*Annulus_Mass.at(1)+Annulus_Mass.at(n_simpson-1);

    for (int j = 2; j < n_simpson; j+=2)
    {
        disk_mass += 2*Annulus_Mass.at(j)+4*Annulus_Mass.at(j+1);
    }
    
    disk_mass *= drad*fourpi/3;
    
    return disk_mass;
}
