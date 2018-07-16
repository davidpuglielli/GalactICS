#include "galaxy_cuda.h"

__host__ __device__  
double SersicPotentialCUDA(double &r)
{
    printf("Error! Can't use analytic Sersic potential on GPU! (Yet)");
//     if (r == 0) return v_bulge_d*v_bulge_d;
//     
//     double u = r/a_bulge_d, notnsersic = 1/n_sersic_d;
//     double un = pow(u, notnsersic);
//     double temp1 = n_sersic_d*(3 - ppp_d);
//     double temp2 = n_sersic_d*(ppp_d - 2);
//     
//     //L1 and L2 are defined in Terzic & Graham (2005). They are used to calculate
//     //the potential of the deprojected Sersic profile.
//     double L1 = rho_0_d*a_bulge_d*a_bulge_d*n_sersic_d*pow(b_n_d,temp2) *
//                 gsl_sf_gamma_inc(-temp2,b_n_d*un);
//     double L2 = rho_0_d*pow(a_bulge_d,3)*n_sersic_d*pow(b_n_d,-temp1)*
//                 (gamma_comp - gsl_sf_gamma_inc(temp1,b_n_d*un));
//     
//     //if (temp1+1 > b_n*un)
//     //    L2 *= gsl_sf_gamma_inc_Q(temp1,b_n*un);
//     //else
//     //    L2 *= (1 - gsl_sf_gamma_inc_Q(temp1,b_n*notnsersic));
//     
//     //cout << "Sersic  " << r << " " << L1 << " " << L2 << " " << fourpi*(L2/r + L1) << endl;
//     //cout << "        " << rho_0*a_bulge*a_bulge*n_sersic*pow(b_n,temp2) << " " 
//     //     << gsl_sf_gamma_inc(-temp2,b_n*un) << " " << -temp2 << " " << b_n*un << endl;
//     
//     return fourpi*(L2/r + L1);
}

__host__ __device__  
double SersicMassCUDA(double &r)
{
//     double u = r/a_bulge_d, notnsersic = 1/n_sersic_d;
//     double un = pow(u, notnsersic);
//     double temp1 = n_sersic_d*(3 - ppp_d);
//     double temp2 = n_sersic_d*(ppp_d - 2);
//         
//     double sersic_mass = fourpi*rho_0_d*pow(a_bulge_d,3)*n_sersic_d*pow(b_n_d,-temp1) *
//                         (gamma_comp - gsl_sf_gamma_inc(temp1,b_n_d*un));
//     
//     return sersic_mass;
}

__host__ __device__  
double SersicDensProfileCUDA(double &r, double *profile)
{
    double u = r / a_bulge_d, notnsersic = 1 / n_sersic_d;
    double un = pow(u, notnsersic);
    
    profile[0] = rho_0_d*pow(u,-ppp_d)*exp(-b_n_d*un);
    
    double deno = a_bulge_d*n_sersic_d*u;
    profile[1] = -SersicDensCUDA(r)*(ppp_d*n_sersic_d + b_n_d*un)/deno;
    
    deno = a_bulge_d*a_bulge_d*n_sersic_d*n_sersic_d*u*u;
    profile[2] = SersicDensCUDA(r)*(ppp_d*ppp_d*n_sersic_d*n_sersic_d + 
                 ppp_d*n_sersic_d*n_sersic_d + 2*ppp_d*b_n_d*n_sersic_d*un + 
                 b_n_d*un*(n_sersic_d - 1) + b_n_d*b_n_d*un*un) / deno;
}

__host__ __device__  
double SersicDensCUDA(double &r)
{
    double u = r/a_bulge_d, notnsersic = 1/n_sersic_d;
    double un = pow(u, notnsersic);
    
    double sersicdens = rho_0_d*pow(u,-ppp_d)*exp(-b_n_d*un);
    
    return sersicdens;
}

__host__ __device__  
double SersicDensPrimeCUDA(double &r)
{
    double u = r / a_bulge_d, notnsersic = 1 / n_sersic_d;
    double un = pow(u, notnsersic);
    
    double deno = a_bulge_d*n_sersic_d*u;
    double sersicdensprime = -SersicDensCUDA(r)*(ppp_d*n_sersic_d + b_n_d*un)/deno;
    
    return sersicdensprime;
}

__host__ __device__  
double SersicDens2PrimeCUDA(double &r)
{
    double u = r / a_bulge_d, notnsersic = 1 / n_sersic_d;
    double un = pow(u, notnsersic);
    
    double deno = a_bulge_d*a_bulge_d*n_sersic_d*n_sersic_d*u*u;
    
    double sersicdens2prime = SersicDensCUDA(r)*
                              (ppp_d*ppp_d*n_sersic_d*n_sersic_d + ppp_d*n_sersic_d*n_sersic_d +
                              2*ppp_d*b_n_d*n_sersic_d*un + b_n_d*un*(n_sersic_d - 1) +
                              b_n_d*b_n_d*un*un) / deno;
    
    return sersicdens2prime;
}
                              
__host__ __device__  
double SersicForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA,
                       thrust::device_ptr<double>B_FR_CUDA)
{
///***  
    if (sersic_flag_d)
    {  
//         double u = r / a_bulge_d, notnsersic = 1 / n_sersic_d;
//         double un = pow(u, notnsersic);
//         double temp1 = n_sersic_d*(3 - ppp_d);
// 
//         double L2 = ppp_d*pow(a_bulge_d,3)*n_sersic_d*pow(b_n_d,-temp1)*
//                     (gamma_comp - gsl_sf_gamma_inc(temp1,b_n_d*un));
// 
//         //if (temp1 > b_n_d*notnsersic)
//         //    L2 *= gsl_sf_gamma_inc(temp1,b_n_d*notnsersic);
//         //else
//         //    L2 *= (1 - gsl_sf_gamma_inc(temp1,b_n_d*notnsersic));
// 
//         return -fourpi*L2/(r*r);
    }
//***/    
   
///*** 
    else   
    {
        double log_r = log10(r);
        //int ihi = (int)((log_r-log_dr)/delta_logr+1);// /(log_rmax-log_dr)*(nr-1)) + 1;
        int ihi = ceil(r/dr_d);

        if(r < dr_d)
        {
            ihi = 1;
        }
        else if (ihi < 1)
        {
        }
        else if (ihi > nr_d-1)
        {
            ihi = nr_d-1;
        }

        double r1 = Radius_CUDA[ihi-1];
        double r2 = Radius_CUDA[ihi];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        return t*B_FR_CUDA[ihi] + tm1*B_FR_CUDA[ihi-1];
    }
//***/    
}
