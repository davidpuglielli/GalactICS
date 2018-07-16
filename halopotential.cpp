#include "galaxy.h"

void HaloPotential(void)
{
    lmax = l_max;
    
    //KK: Now get the harmonics of the density in this potential --> Halo_Dens
    int n_theta = lmax*4 + 2;
    double eps = 0.0001;
    Halo_Dens[0][0] = HaloDens(eps, 0) * sqrt4pi;
    
    for (int l = 2; l < lmax+1; l+=2)
        Halo_Dens[l/2][0] = 0;
    
    for (int ir = 0; ir < nr; ++ir)
        Halo_Dens[0][ir] = 0;
    
    //KK: nrmx will mark the outermost radial bin with non-zero density.
    int nrmx = nr;
    
    for (int l = 0; l < lmax+1; l+=2)
    {
        //KK: integrate density * spherical harmonic function over quadrant
        for (int ir = 0; ir < nrmx; ++ir)
        {//cout<<"heeere"<<endl;
            double rad = Radius[ir];
            double s = 0;
            double d_cos_theta = 1.0/n_theta;
            
            s += PolarHaloDens(rad, 1, l) + PolarHaloDens(rad, 0, l);
            
            for (int is = 1; is < n_theta; is+=2)
            {
                double cos_theta = is*d_cos_theta;
                s += 4*PolarHaloDens(rad, cos_theta, l);
            }
            
            for (int is = 2; is < n_theta-1; is+=2)
            {
                double cos_theta = is*d_cos_theta;
                s += 2*PolarHaloDens(rad, cos_theta, l);
            }
            
            s *= d_cos_theta*4*PI/3;
            Halo_Dens[l/2][ir] = s;
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
            S_1[ir] = S_1[ir-2] + dr/3*(Halo_Dens[l/2][ir-2]*pow(rad-2*dr, l+2) +
                     4*Halo_Dens[l/2][ir-1]*pow(rad-dr,l+2) + Halo_Dens[l/2][ir]*pow(rad,l+2));
        }
        
        S_2[nr-1] = 0;
        
        for (int ir = nr-3; ir > 1; ir-=2)
        {
            double rad = Radius[ir];
            S_2[ir] = S_2[ir+2] + dr/3*(Halo_Dens[l/2][ir+2]*pow(rad+2*dr, 1-l) +
                     4*Halo_Dens[l/2][ir+1]*pow(rad+dr,1-l) + Halo_Dens[l/2][ir]*pow(rad,1-l));
        }
        
        for (int ir = 2; ir < nr; ir+=2)
        {
            double rad = Radius[ir];
            Halo_Pot[l/2][ir] = -fourpi * (S_1[ir]/pow(rad,l+1) + S_2[ir]*pow(rad,l)) / (2*l+1);
        }
        
    //KK: Calculate the radial gradients
        for (int ir = 2; ir < nr; ir+=2)
        {
            double rad = Radius[ir];
            Halo_FR[l/2][ir] = -fourpi*(-(l+1)*S_1[ir]/pow(rad,l+2) +
                              l*S_2[ir]*pow(rad,l-1)) / (2.0*l+1);
        }
    }
        
     //KK: now interpolate the gaps first quadratically interpolate the
     //monopole back to the origin.  the remaining multipoles are zero
     //there.
        
    Halo_Pot[0][0] = 3*(A_Pot[0][2] - A_Pot[0][4]) + A_Pot[0][6];
    Halo_FR[0][0] = 0;

    for (int l = 2; l < lmax+1; l+=2)
    {
        Halo_Pot[l/2][0] = 0;
        Halo_FR[l/2][0] = 0;
    }

    for (int ir = 1; ir < nr; ir+=2)
    {
        for (int l = 0; l < lmax+1; l+=2)
        {
            Halo_Pot[l/2][ir] = 0.5*Halo_Pot[l/2][ir-1] + 0.5*Halo_Pot[l/2][ir+1];
            Halo_FR[l/2][ir] = 0.5*Halo_FR[l/2][ir-1] + 0.5*Halo_FR[l/2][ir+1];
            //cout << "halofr   " << ir << " " << Halo_FR[l/2][ir] << endl;
        }
    }

    double a00 = Halo_Pot[0][0];

    for (int ir = 0; ir < nr; ++ir)
    {
        Halo_Pot[0][ir] -= a00;
    }
                 
 //KK: write halo final potential.
 //Reset the potentials so that phi is 0 at infinity
    
    double r_edge = Radius[nr-1];
    
    for (int l = 0; l < lmax+1; l+=2)
    {
        double constant;//cout << "lmax " << lmax << " " <<Halo_Pot[l/2][nr-1]<<" "<<Halo_FR[l/2][nr-1]<<endl; 
        constant = Halo_Pot[l/2][nr-1] + Halo_FR[l/2][nr-1]*r_edge/(l+1);
        
        for (int i = 0; i < nr; ++i)
        {
            Halo_Pot[l/2][i] = Halo_Pot[l/2][i] - constant;
        }
    }
    
    if (do_file_io)
    {
        WriteHDat();
    }
    
    halo_mass = Halo_FR[0][nr-1]*Radius[nr-1]*Radius[nr-1]*oneoversqrt4pi;
    
    for (int i = 0; i < nr; ++i)
    {
        if(A_Dens[0][i] == 0)
        {
            halo_edge=Radius[i];
            break;
        }
    }
    
    cout << "Halo mass = " << halo_mass << endl;
    cout << "Halo edge radius = " << halo_edge << endl;
}
    
