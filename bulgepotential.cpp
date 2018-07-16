#include "galaxy.h"

void BulgePotential(void)
{
    lmax = l_max;
    
    //KK: Now get the harmonics of the density in this potential --> Bulge_Dens
    int n_theta = lmax*4 + 2;
    double eps = 0.00001;
    Bulge_Dens[0][0] = BulgeDens(eps, 0) * sqrt4pi;
    
    for (int l = 2; l < lmax+1; l+=2)
        Bulge_Dens[l/2][0] = 0;
    
    for (int ir = 0; ir < nr; ++ir)
        Bulge_Dens[0][ir] = 0;
    
    //KK: nrmx will mark the outermost radial bin with non-zero density.
    int nrmx = nr;
    
    for (int l = 0; l < lmax+1; l+=2)
    {
        for (int ir = 0; ir < nrmx; ++ir)
        {
            double rad = Radius[ir];
            double s = 0;
            double d_cos_theta = 1.0/n_theta;
            
            s += PolarBulgeDens(rad, 1, l) + PolarBulgeDens(rad, 0, l);
            
            for (int is = 1; is < n_theta; is+=2)
            {
                double cos_theta = is*d_cos_theta;
                s += 4*PolarBulgeDens(rad, cos_theta, l);
            }
            
            for (int is = 2; is < n_theta-1; is+=2)
            {
                double cos_theta = is*d_cos_theta;
                s += 2*PolarBulgeDens(rad, cos_theta, l);
            }
            
            s *= d_cos_theta*fourpi*0.333333333333;
            Bulge_Dens[l/2][ir] = s;
        }
    }
        
    //KK: now get the potential harmonics of this new density. Simpson's
    //rule integration. (BT 2-208)
    //DP: That would be the *first* edition of BT
        
    for (int l = 0; l < lmax+1; l+=2)
    {
        S_1[0] = 0;
        
        for (int ir = 2; ir < nr; ir+=2)
        {
            double rad = Radius[ir];
            S_1[ir] = S_1[ir-2] + dr/3*(Bulge_Dens[l/2][ir-2]*pow(rad-2*dr,l+2) +
                     4*Bulge_Dens[l/2][ir-1]*pow(rad-dr,l+2) + Bulge_Dens[l/2][ir]*pow(rad,l+2));
            //cout << "s1  " << rad << " " << S_1[ir] << endl;
        }
        
        S_2[nr-1] = 0;
        
        for (int ir = nr-3; ir > 1; ir-=2)
        {
            double rad = Radius[ir];
            S_2[ir] = S_2[ir+2] + dr/3*(Bulge_Dens[l/2][ir+2]*pow(rad+2*dr, 1-l) +
                     4*Bulge_Dens[l/2][ir+1]*pow(rad+dr,1-l) + Bulge_Dens[l/2][ir]*pow(rad,1-l));
            //cout << "s2  " << rad << " " << S_2[ir] << endl;
        }
        
        for (int ir = 2; ir < nr; ir+=2)
        {
            double rad = Radius[ir];
            Bulge_Pot[l/2][ir] = -fourpi * (S_1[ir]/pow(rad,l+1) + S_2[ir]*pow(rad,l)) / (2*l+1);
        }
        
    //KK: Calculate the radial gradients
        for (int ir = 2; ir < nr; ir+=2)
        {
            double rad = Radius[ir];
            Bulge_FR[l/2][ir] = -fourpi*(-(l+1)*S_1[ir]/pow(rad,l+2) +
                                l*S_2[ir]*pow(rad,l-1)) / (2.0*l+1);
        }
    }
        
     //KK: now interpolate the gaps first quadratically interpolate the
     //monopole back to the origin.  the remaining multipoles are zero
     //there.
        
    Bulge_Pot[0][0] = 3*(Bulge_Pot[0][2] - Bulge_Pot[0][4]) + Bulge_Pot[0][6];
    Bulge_FR[0][0] = 0;

    for (int l = 2; l < lmax+1; l+=2)
    {
        Bulge_Pot[l/2][0] = 0;
        Bulge_FR[l/2][0] = 0;
    }

    for (int ir = 1; ir < nr; ir+=2)
    {
        for (int l = 0; l < lmax+1; l+=2)
        {
            Bulge_Pot[l/2][ir] = 0.5*Bulge_Pot[l/2][ir-1] + 0.5*Bulge_Pot[l/2][ir+1];
            Bulge_FR[l/2][ir] = 0.5*Bulge_FR[l/2][ir-1] + 0.5*Bulge_FR[l/2][ir+1];
            //cout << "potfr  " << Radius[ir] << " " << Bulge_Pot[l/2][ir] << " " 
            //     << Bulge_FR[l/2][ir] << endl;
        }
    }

    double a00 = Bulge_Pot[0][0];

    for (int ir = 0; ir < nr; ++ir)
    {
        Bulge_Pot[0][ir] -= a00;
    }
                 
 //KK: write bulge final potential.
 //Reset the potentials so that phi is 0 at infinity
    
    r_edge = Radius[nr-1];
    
    for (int l = 0; l < lmax+1; l+=2)
    {
        double constant;
        constant = Bulge_Pot[l/2][nr-1] + Bulge_FR[l/2][nr-1]*r_edge/(l+1);
        
        for (int i = 0; i < nr; ++i)
        {
            Bulge_Pot[l/2][i] = Bulge_Pot[l/2][i] - constant;
        }
    }
    
    //write to b.dat
    if (do_file_io)
    {
        WriteBDat();
    }
    
    bulge_mass = Bulge_FR[0][nr-1]*Radius[nr-1]*Radius[nr-1]*oneoversqrt4pi;
    
    for (int i = 0; i < nr; ++i)
    {
        if(A_Dens[0][i] == 0)
        {
            bulge_edge=Radius[i];
            break;
        }
    }
    
    cout << "Bulge mass = " << bulge_mass << endl;
    cout << "Bulge edge radius = " << bulge_edge << endl;
}
    
void BulgePotentialEstimate(void)
{
    double eps = 0.0001;
    
    for (int i = 0; i < nr; ++i)
    {
        double r = Radius[i];
        if (r == 0)
            r = eps;
        B_Dens[i]=SersicDens(r);
#ifdef DEBUG
        //cout << i << " " << r << " " << HaloDensity(r) << endl;
#endif
    }
    
// KK: now get the potential harmonics of this new density. (BT 2-208)
// Simpson's rule integration.

    double ds = delta_logr;
    
    S_1[0] = 0;
    double r = Radius[2];
    double ss = log10(r);
    
    S_1[2] = (r*dr/3)*(4*B_Dens[1]*(1-dr/r)*(1-dr/r)+B_Dens[2]);
    //S_1[2] = (ds/3) * (4*B_Dens[1]*exp(3*ss-3*ds) + B_Dens[2]*exp(3*ss));cout<<"s1   "<<S_1[2]<<endl;
    
    double r_old = r;
    
    for (int ir = 4; ir < nr; ir+=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s1a = (ds/3)*(B_Dens[ir-2]*exp(3*s-6*ds) + 4*B_Dens[ir-1]*exp(3*s-3*ds) +
        //                     B_Dens[ir]*exp(3*s));
        double s1a = (r*dr/3)*(B_Dens[ir-2]*pow(1-2*dr/r, 2) + 
                     4*B_Dens[ir-1]*pow(1-dr/r, 2) + B_Dens[ir]);
        S_1[ir] = s1a+S_1[ir-2]*r_old/r;
        r_old = r;
        //cout<<"s1   "<<r<<" "<<s1a<<" "<<S_1[ir]<<endl;
    }
    
    S_2[nr-1] = 0;
    
    r_old = Radius[nr-1];
    
    for (int ir = nr-3; ir > 1; ir-=2)
    {
        r = Radius[ir];
        double s = log10(r);
        
        //double s2a = (ds/3)*(B_Dens[ir]*exp(3*s) + 4*B_Dens[ir+1]*exp(3*s+3*ds) +
        //                     B_Dens[ir+2]*exp(3*s+6*ds));
        double s2a = (r*dr/3) * (B_Dens[ir+2]*(1+2*dr/r) + 
                     4*B_Dens[ir+1]*(1+dr/r) + B_Dens[ir]);
        S_2[ir] = s2a + S_2[ir+2];
        //cout<<"s2 h   "<<r<<" "<<s2a<<" "<<S_2[ir]<<endl;
        r_old = r;
    }
    
    for (int ir = 2; ir < nr; ir+=2)
    {
        r = Radius[ir];
        B_Pot[ir] = fourpi*(S_1[ir]+S_2[ir]);
        B_FR[ir] = -fourpi*S_1[ir]/r;
#ifdef DEBUG
        //cout << "Bulgepotentialestimate  " << r << " " << ir << " " << B_Pot[ir] << " " << B_FR[ir] << endl;
#endif
    }
    
    B_Pot[0] = 3*(B_Pot[2]-B_Pot[4])+B_Pot[6];
    //cout << "Halopotentialestimate  " << r << " " << B_Pot[0] << " " << B_FR[0] << endl;
    B_FR[0] = 0;
    
//  KK: then linearly interpolate other bins.
    
    for (int ir = 1; ir < nr; ir+=2)
    {
        B_Pot[ir] = (B_Pot[ir-1]+B_Pot[ir+1])/2;
        B_FR[ir] = (B_FR[ir-1]+B_FR[ir+1])/2;
#ifdef DEBUG
        //cout << "Bulgepotentialestimate  " << Radius[ir] << " " << ir << " " << B_Pot[ir] << " " << B_FR[ir] << endl;
#endif
    }
}
