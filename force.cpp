#include "galaxy.h"

void Force(double &s, double &z, double &fs, double &fz)
{
    //double param[5] = {0.282094792, 0.630783131, 0.846284375, 
    //                   1.017107236, 1.163106623};
    double pc[20], p[20], dp[20], pot;
    
    double r = sqrt(s*s + z*z);
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    int ihim = ihi - 1;
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Force finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<s<<" "<<z<<" "<< endl;
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
    double r_edge = Radius[nr-1];
    double t = (r-r1)/(r2-r1);
    double tm = 1 - t;
    
    if (r == 0)
    {
        fs = fz = 0;
    }
    else
    {
        double cos_theta = z/r;
        double cos2theta = cos_theta*cos_theta;
        double sin_theta = s/r;
        double sin2theta = 1 - cos2theta;
        
        for (int l = 0; l < lmax+1; l+=2)
        {
            pc[l/2] = sqrt((2.0*l+1)/4/PI);
            p[l/2] = gsl_sf_legendre_Pl(l, cos_theta);
            
            if (fabs(cos_theta) == 1)
            {
                dp[l/2] = 0;
            }
            else
            {
                dp[l/2] = l*(gsl_sf_legendre_Pl(l, cos_theta)-cos_theta*p[l/2]) / sin2theta;
                //cout<<"dp "<<dp[l/2]<<" "<<sin2theta<<" "<<cos_theta<<endl;
            }   
        }
        
        for (int i = 0; i < lmax/2 + 1; ++i)
        {
            p[i] = p[i]*pc[i];
            dp[i] = dp[i]*pc[i];
        }
        
        double frr = 0, fth = 0, pot = 0;
        
        if (r <= r_edge)
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                frr += p[i]*(t*F_R[i][ihi] + tm*F_R[i][ihim]);
                //cout << "force1 " << frr << " " << fth << " " << lmax << endl;
            }
            //cout<<"here "<<r<<" "<<ihi<<" "<<ihim<<" "<<dp[2-1]<<endl;
            //cout<<"here "<<" "<<sin_theta<<" "<<A_Pot[2][ihim]
            //        <<" "<<A_Pot[2][ihim]<<" "<<dp[2-1]<<endl;
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                //cout<<"here "<<i<<" "<<sin_theta<<" "<<A_Pot[i][ihim]
                //    <<" "<<A_Pot[i-1][ihim]<<" "<<dp[i-1]<<endl;
                fth -= sin_theta*dp[i-1]*(t*A_Pot[i-1][ihi] + tm*A_Pot[i-1][ihim]);
                //cout<<"here "<<i<<" " << r<<" " <<z<<" "<<A_Pot[i-1][ihi]<<" "
                //    <<A_Pot[i-1][ihim]<<"   "<<t<<" "<<tm<<" "<<sin_theta<<" "
                //    <<dp[i-1]<<" "<<fth<< endl;
            }
             //cout<<"here "<<endl;
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                pot += p[i]*(t*A_Pot[i][ihi] + tm*A_Pot[i][ihim]);
            }
            
            //cout << "force " << frr << " " << fth << endl;
        }
        else
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                frr -= (l+1)*p[i]*A_Pot[i][nr] / r_edge * pow(r_edge/r,l+2);
                //cout << "force2 " << r << " " << r_edge << " " << frr << " " << fth << " " << lmax << endl;
            }
            
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                int l = 2*(i-1);
                fth -= sin_theta*dp[i-1]*A_Pot[i-1][nr] * pow(r_edge/r,l+1);
            }
            
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                pot += p[i]*A_Pot[i][nr] * pow(r_edge/r,l+1);
            }
        }
        
        if (disk_flag == 1)
        {
            pot += AppDiskPot(s,z);
        }
        
        fs = -sin_theta*frr - cos_theta*fth/r;
        fz = -cos_theta*frr + sin_theta*fth/r;
        //cout << "force " << frr << " " << fth << " " << lmax << endl;        
        if (disk_flag == 1)
        {
            double fsad, fzad;
            AppDiskForce(s,z,fsad,fzad);
            fs += fsad;
            fz += fzad;
        }
    }
    
    //cout << "force " << fs << " " << fz << endl;
}    

void HaloForce(double &s, double &z, double &fs, double &fz)
{
    //double param[5] = {0.282094792, 0.630783131, 0.846284375, 
    //                   1.017107236, 1.163106623};
    double pc[20], p[20], dp[20], pot;
    
    double r = sqrt(s*s + z*z);
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    int ihim = ihi - 1;
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Force finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<s<<" "<<z<<" "<< endl;
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
    double r_edge = Radius[nr-1];
    double t = (r-r1)/(r2-r1);
    double tm = 1 - t;
    
    if (r == 0)
    {
        fs = fz = 0;
    }
    else
    {
        double cos_theta = z/r;
        double cos2theta = cos_theta*cos_theta;
        double sin_theta = s/r;
        double sin2theta = 1 - cos2theta;
        
        for (int l = 0; l < lmax+1; l+=2)
        {
            pc[l/2] = sqrt((2.0*l+1)/4/PI);
            p[l/2] = gsl_sf_legendre_Pl(l, cos_theta);
            
            if (fabs(cos_theta) == 1)
            {
                dp[l/2] = 0;
            }
            else
            {
                dp[l/2] = l*(gsl_sf_legendre_Pl(l, cos_theta)-cos_theta*p[l/2]) / sin2theta;
                //cout<<"dp "<<dp[l/2]<<" "<<sin2theta<<" "<<cos_theta<<endl;
            }   
        }
        
        for (int i = 0; i < lmax/2 + 1; ++i)
        {
            p[i] = p[i]*pc[i];
            dp[i] = dp[i]*pc[i];
        }
        
        double frr = 0, fth = 0, pot = 0;
        
        if (r <= r_edge)
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                frr += p[i]*(t*Halo_FR[i][ihi] + tm*Halo_FR[i][ihim]);
                //cout << "force1 " << frr << " " << fth << " " << lmax << endl;
            }
            //cout<<"here "<<r<<" "<<ihi<<" "<<ihim<<" "<<dp[2-1]<<endl;
            //cout<<"here "<<" "<<sin_theta<<" "<<A_Pot[2][ihim]
            //        <<" "<<A_Pot[2][ihim]<<" "<<dp[2-1]<<endl;
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                //cout<<"here "<<i<<" "<<sin_theta<<" "<<A_Pot[i][ihim]
                //    <<" "<<A_Pot[i-1][ihim]<<" "<<dp[i-1]<<endl;
                fth -= sin_theta*dp[i-1]*(t*Halo_Pot[i-1][ihi] + tm*Halo_Pot[i-1][ihim]);
                //cout<<"here "<<i<<" " << r<<" " <<z<<" "<<A_Pot[i-1][ihi]<<" "
                //    <<A_Pot[i-1][ihim]<<"   "<<t<<" "<<tm<<" "<<sin_theta<<" "
                //    <<dp[i-1]<<" "<<fth<< endl;
            }
             //cout<<"here "<<endl;
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                pot += p[i]*(t*Halo_Pot[i][ihi] + tm*Halo_Pot[i][ihim]);
            }
            
            //cout << "force " << frr << " " << fth << endl;
        }
        else
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                frr -= (l+1)*p[i]*Halo_Pot[i][nr] / r_edge * pow(r_edge/r,l+2);
                //cout << "force2 " << r << " " << r_edge << " " << frr << " " << fth << " " << lmax << endl;
            }
            
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                int l = 2*(i-1);
                fth -= sin_theta*dp[i-1]*Halo_Pot[i-1][nr] * pow(r_edge/r,l+1);
            }
            
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                pot += p[i]*Halo_Pot[i][nr] * pow(r_edge/r,l+1);
            }
        }
        
        if (disk_flag == 1)
        {
            pot += AppDiskPot(s,z);
        }
        
        fs = -sin_theta*frr - cos_theta*fth/r;
        fz = -cos_theta*frr + sin_theta*fth/r;
        //cout << "force " << frr << " " << fth << " " << lmax << endl;        
        if (disk_flag)
        {
            double fsad, fzad;
            AppDiskForce(s,z,fsad,fzad);
            fs += fsad;
            fz += fzad;
        }
    }
    
    //cout << "force " << fs << " " << fz << endl;
}   
 
void BulgeForce(double &s, double &z, double &fs, double &fz)
{
    //double param[5] = {0.282094792, 0.630783131, 0.846284375, 
    //                   1.017107236, 1.163106623};
    double pc[20], p[20], dp[20], pot;
    
    double r = sqrt(s*s + z*z);
    double log_r = log10(r);
    //int ihi = (int)((log_r-log_dr)/(log_rmax-log_dr)*(nr-1)) + 1;
    int ihi = ceil(r/dr);
    int ihim = ihi - 1;
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Force finds out of range indices. Exiting..." << endl;
        cout <<"ihi "<<ihi<<" "<<r<<" "<<s<<" "<<z<<" "<< endl;
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
    double r_edge = Radius[nr-1];
    double t = (r-r1)/(r2-r1);
    double tm = 1 - t;
    
    if (r == 0)
    {
        fs = fz = 0;
    }
    else
    {
        double cos_theta = z/r;
        double cos2theta = cos_theta*cos_theta;
        double sin_theta = s/r;
        double sin2theta = 1 - cos2theta;
        
        for (int l = 0; l < lmax+1; l+=2)
        {
            pc[l/2] = sqrt((2.0*l+1)/4/PI);
            p[l/2] = gsl_sf_legendre_Pl(l, cos_theta);
            
            if (fabs(cos_theta) == 1)
            {
                dp[l/2] = 0;
            }
            else
            {
                dp[l/2] = l*(gsl_sf_legendre_Pl(l, cos_theta)-cos_theta*p[l/2]) / sin2theta;
                //cout<<"dp "<<dp[l/2]<<" "<<sin2theta<<" "<<cos_theta<<endl;
            }   
        }
        
        for (int i = 0; i < lmax/2 + 1; ++i)
        {
            p[i] = p[i]*pc[i];
            dp[i] = dp[i]*pc[i];
        }
        
        double frr = 0, fth = 0, pot = 0;
        
        if (r <= r_edge)
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                frr += p[i]*(t*Bulge_FR[i][ihi] + tm*Bulge_FR[i][ihim]);
                //cout << "force1 " << frr << " " << fth << " " << lmax << endl;
            }
            //cout<<"here "<<r<<" "<<ihi<<" "<<ihim<<" "<<dp[2-1]<<endl;
            //cout<<"here "<<" "<<sin_theta<<" "<<A_Pot[2][ihim]
            //        <<" "<<A_Pot[2][ihim]<<" "<<dp[2-1]<<endl;
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                //cout<<"here "<<i<<" "<<sin_theta<<" "<<A_Pot[i][ihim]
                //    <<" "<<A_Pot[i-1][ihim]<<" "<<dp[i-1]<<endl;
                fth -= sin_theta*dp[i-1]*(t*Bulge_Pot[i-1][ihi] + tm*Bulge_Pot[i-1][ihim]);
                //cout<<"here "<<i<<" " << r<<" " <<z<<" "<<A_Pot[i-1][ihi]<<" "
                //    <<A_Pot[i-1][ihim]<<"   "<<t<<" "<<tm<<" "<<sin_theta<<" "
                //    <<dp[i-1]<<" "<<fth<< endl;
            }
             //cout<<"here "<<endl;
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                pot += p[i]*(t*Bulge_Pot[i][ihi] + tm*Bulge_Pot[i][ihim]);
            }
            
            //cout << "force " << frr << " " << fth << endl;
        }
        else
        {
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                frr -= (l+1)*p[i]*Bulge_Pot[i][nr] / r_edge * pow(r_edge/r,l+2);
                //cout << "force2 " << r << " " << r_edge << " " << frr << " " << fth << " " << lmax << endl;
            }
            
            for (int i = 2; i < lmax/2 + 2; ++i)
            {
                int l = 2*(i-1);
                fth -= sin_theta*dp[i-1]*Bulge_Pot[i-1][nr] * pow(r_edge/r,l+1);
            }
            
            for (int i = 0; i < lmax/2 + 1; ++i)
            {
                int l = 2*i;
                pot += p[i]*Bulge_Pot[i][nr] * pow(r_edge/r,l+1);
            }
        }
        
        if (disk_flag == 1)
        {
            pot += AppDiskPot(s,z);
        }
        
        fs = -sin_theta*frr - cos_theta*fth/r;
        fz = -cos_theta*frr + sin_theta*fth/r;
        //cout << "force " << frr << " " << fth << " " << lmax << endl;        
        if (disk_flag == 1)
        {
            double fsad, fzad;
            AppDiskForce(s,z,fsad,fzad);
            fs += fsad;
            fz += fzad;
        }
    }
    
    //cout << "force " << fs << " " << fz << endl;
}    
