//These functions calculate quantities related to sersic profiles. The relevant
//expressions are found in Prugniel & Simien (1997) and Terzic & Graham (2005)
//and GSL routines are used to compute the various Gamma functions.

#include "galaxy.h"

double SersicPotential(double &r)
{
    if (r == 0) return G.v_bulge*G.v_bulge;
    
    double u = r/G.a_bulge, notnsersic = 1/G.n_sersic;
    double un = pow(u, notnsersic);
    double temp1 = G.n_sersic*(3 - G.ppp);
    double temp2 = G.n_sersic*(G.ppp - 2);
    
    //L1 and L2 are defined in Terzic & Graham (2005). They are used to calculate
    //the potential of the deprojected Sersic profile.
    double L1 = G.rho_0*G.a_bulge*G.a_bulge*G.n_sersic*pow(G.b_n,temp2) *
                gsl_sf_gamma_inc(-temp2,G.b_n*un);
    double L2 = G.rho_0*pow(G.a_bulge,3)*G.n_sersic*pow(G.b_n,-temp1)*
                (gamma_comp - gsl_sf_gamma_inc(temp1,G.b_n*un));
    
    //if (temp1+1 > G.b_n*un)
    //    L2 *= gsl_sf_gamma_inc_Q(temp1,G.b_n*un);
    //else
    //    L2 *= (1 - gsl_sf_gamma_inc_Q(temp1,G.b_n*notnsersic));
    
    //cout << "Sersic  " << r << " " << L1 << " " << L2 << " " << fourpi*(L2/r + L1) << endl;
    //cout << "        " << G.rho_0*G.a_bulge*G.a_bulge*G.n_sersic*pow(G.b_n,temp2) << " " 
    //     << gsl_sf_gamma_inc(-temp2,G.b_n*un) << " " << -temp2 << " " << G.b_n*un << endl;
    
    return fourpi*(L2/r + L1);
}

double SersicMass(double &r)
{
    double u = r/G.a_bulge, notnsersic = 1/G.n_sersic;
    double un = pow(u, notnsersic);
    double temp1 = G.n_sersic*(3 - G.ppp);
    double temp2 = G.n_sersic*(G.ppp - 2);
        
    double sersic_mass = fourpi*G.rho_0*pow(G.a_bulge,3)*G.n_sersic*pow(G.b_n,-temp1) *
                        (gamma_comp - gsl_sf_gamma_inc(temp1,G.b_n*un));
    
    return sersic_mass;
}

double SersicDensProfile(double &r, double *profile)
{
    double u = r / G.a_bulge, notnsersic = 1 / G.n_sersic;
    double un = pow(u, notnsersic);
    
    profile[0] = G.rho_0*pow(u,-G.ppp)*exp(-G.b_n*un);
    
    double deno = G.a_bulge*G.n_sersic*u;
    profile[1] = -SersicDens(r)*(G.ppp*G.n_sersic + G.b_n*un)/deno;
    
    deno = G.a_bulge*G.a_bulge*G.n_sersic*G.n_sersic*u*u;
    profile[2] = SersicDens(r)*(G.ppp*G.ppp*G.n_sersic*G.n_sersic + 
                 G.ppp*G.n_sersic*G.n_sersic + 2*G.ppp*G.b_n*G.n_sersic*un + 
                 G.b_n*un*(G.n_sersic - 1) + G.b_n*G.b_n*un*un) / deno;
}

double SersicDens(double &r)
{
    double u = r/G.a_bulge, notnsersic = 1/G.n_sersic;
    double un = pow(u, notnsersic);
    
    double sersicdens = G.rho_0*pow(u,-G.ppp)*exp(-G.b_n*un);
    
    return sersicdens;
}

double SersicDensPrime(double &r)
{
    double u = r / G.a_bulge, notnsersic = 1 / G.n_sersic;
    double un = pow(u, notnsersic);
    
    double deno = G.a_bulge*G.n_sersic*u;
    double sersicdensprime = -SersicDens(r)*(G.ppp*G.n_sersic + G.b_n*un)/deno;
    
    return sersicdensprime;
}

double SersicDens2Prime(double &r)
{
    double u = r / G.a_bulge, notnsersic = 1 / G.n_sersic;
    double un = pow(u, notnsersic);
    
    double deno = G.a_bulge*G.a_bulge*G.n_sersic*G.n_sersic*u*u;
    
    double sersicdens2prime = SersicDens(r)*
                              (G.ppp*G.ppp*G.n_sersic*G.n_sersic + G.ppp*G.n_sersic*G.n_sersic +
                              2*G.ppp*G.b_n*G.n_sersic*un + G.b_n*un*(G.n_sersic - 1) +
                              G.b_n*G.b_n*un*un) / deno;
    
    return sersicdens2prime;
}
                              
double SersicForce(double &r)
{
///***  
    if (sersic_flag)
    {  
        double u = r / G.a_bulge, notnsersic = 1 / G.n_sersic;
        double un = pow(u, notnsersic);
        double temp1 = G.n_sersic*(3 - G.ppp);

        double L2 = G.rho_0*pow(G.a_bulge,3)*G.n_sersic*pow(G.b_n,-temp1)*
                    (gamma_comp - gsl_sf_gamma_inc(temp1,G.b_n*un));

        //if (temp1 > G.b_n*notnsersic)
        //    L2 *= gsl_sf_gamma_inc(temp1,G.b_n*notnsersic);
        //else
        //    L2 *= (1 - gsl_sf_gamma_inc(temp1,G.b_n*notnsersic));

        return -fourpi*L2/(r*r);
    }
//***/    
   
///*** 
    else   
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
            cout << "SersicForce finds out of range indices. Exiting..." << endl;
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

        return t*B_FR[ihi] + tm1*B_FR[ihi-1];
    }
//***/    
}

//getting the coefficient b which is dependent on the Sersic index n_sersic.
//This is the approximation found by Prugniel & simien. The more exact solution
//which involves solving an equation with Gamma functions will be implemented later
double Get_b_n(void)
{
    return 2*G.n_sersic-0.33333333+0.009876/G.n_sersic;
}

//The gamma function in this function needs to be fixed
double Get_rho_0(void)
{
    double temp = G.n_sersic*(G.ppp - 2), fixthis=1;
    double den = G.a_bulge*G.a_bulge*G.n_sersic*pow(G.b_n,temp)*gsl_sf_gamma(-temp);
    
    return G.v_bulge*G.v_bulge*oneover4pi/den;
}
      
    
    
