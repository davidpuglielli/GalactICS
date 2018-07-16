#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "galaxy.h"

//put the global variables here? Update them on the host?
__constant__ int n_psi_d, n_int_d, nr_d, nr_spline_d, n_iter_d, nrmx_d, lmax_d, l_max_d; 
__constant__ int max_iter_d, iter_d, n_simpson_d;
__constant__ int halo_flag_d, disk_flag_d, gasdisk_flag_d, bulge_flag_d, smbh_flag_d;
__constant__ int do_file_io_d, nbody_flag_d, chisq_flag_d, sersic_flag_d;
__constant__ int disk_d, gas_disk_d, gasdisk_params_d, disk_params_d, nondisk_params_d;
__constant__ int astro_params_d, error_params_d;

__constant__ double potcor_d, psi_crit_d, psi_d_d, psi_0_d, tolerance_factor_d, fraction_d;
__constant__ double psi_upper_initial_d, psi_lower_initial_d, gamma_comp_d, log_rj2_d;
__constant__ double dr_d, r_max_d, log_rmax_d, log_dr_d, delta_logr_d;
__constant__ double total_mass_d, halo_mass_d, disk_mass_d, total_disk_mass_d, bulge_mass_d;
__constant__ double r_edge_d, halo_edge_d, disk_edge_d, bulge_edge_d;
__constant__ double Q_total_d, X_total_d;
__constant__ double sin_inclination_d, cos_inclination_d;

__constant__ double Plcon_d[40];

struct TotalDensityFunctor
{
    const galaxy G;
    //const float b;
    const thrust::device_ptr<float> A_Pot;
    const thrust::device_ptr<float> A_Density;
    
    saxpy_functor(galstruct _G, thrust::device_ptr<float> _A_Pot, 
                                thrust::device_ptr<float> _A_Dens):
            G(_G),A_Pot(_A_Pot),A_Dens(_A_Dens){}
    //saxpy_functor(galstruct _a, float _b):a(_a), b(_b){}
    //saxpy_functor(){}
            
    __host__ __device__ 
    float operator()(const float& rad, const float& y) const
    {
        int n_theta = 100;
        double cos_theta, d_cos_theta = 1.0/n_theta;
        
        for (int j = 0; j < n_theta+1; ++j)
        {
            double costheta = j*d_cos_theta;
            double z = rad*costheta;
            double ss = rad*sqrt(1-cos_theta*cos_theta);

            double psi = PotCUDA(rad,z,A_Pot);

            double totdens = 0;

            for (int i = 0; i < disk; ++i)
            {
                if(fabs(z/G.Z_Disk[i])>30)
                    continue;

                double trunc_fac = GetTrunc(r, i);

                if (trunc_fac==0)
                    continue;

                double con;

                if (z==0)
                {
                    con = 1;
                }
                else if (!halo_flag && !bulge_flag)
                {
                    con = exp(-fabs(z/G.Z_Disk[i]));
                    con = 2.0*con/(1.0+con*con);
                    con *= con;
                }
                else
                {
                    //KK: make disk vertically isothermal; normalize density so that
                    //at height 3*zdisk the vertical density drop is that of a sech^2.
                    double psizh = PotCUDA(r, 3.0*G.Z_Disk[i]);
                    double psi00 = PotCUDA(r, 0);
                    double dpsizh = psi00 - psizh;
                    double dpsi = psi00 - psi;
                    double coeff = 0;

                    if (r==0 && z==0)
                    {
                        dpsi = 0;
                    }

                    if (fabs(dpsizh) > 0)
                    {
                        coeff = dpsi / dpsizh;
                    }

                    if (coeff > 16 || coeff < 0)
                    {
                        con = 0;
                    }
                    else
                    {
                        con = pow(0.009866, coeff);
                    }
                }

	            double Sigma_Profile[2], Rho_Profile[2];

                DiskProfile(r, z, i, Sigma_Profile, Rho_Profile);
                totdens += Rho_Disk_Const[i]*Sigma_Profile[0]*con*trunc_fac;
            }

            //Do the halo            
            if (psi < psi_crit)
            {
                totdens+=0;
            }                    
            else if (psi >= psi_0) 
            {
                totdens += Dens_Psi_Halo[0];
            }
            else
            {
                double rj1 = (psi_0 - psi)/(psi_0 - psi_d);
                //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
                double rj = (n_psi-1.0)*log(rj1)/log(rj2);

                int j = int(rj);

                if (j < 0)
                    j = 0;
                else if (j >= n_psi-1)
                    j = n_psi - 2;

                double frac = rj - j;

                double halo_dens_psi = *(Dens_Psi_Halo_CUDA+j) + //Dens_Psi_Halo[j] + 
                                     frac*(*(Dens_Psi_Halo_CUDA+j+1)-(*(Dens_Psi_Halo_CUDA+j)));
                                     //  frac*(Dens_Psi_Halo[j+1]-Dens_Psi_Halo[j]);
                totdens += halo_dens_psi;
            }

            //Do the bulge
            if (psi < psi_crit)
            {
                totdens+=0;
            }                    
            else if (psi >= psi_0) 
            {
                totdens += *(Dens_Psi_Bulge_CUDA);
            }
            else
            {
                double rj1 = (psi_0 - psi)/(psi_0 - psi_d);
                //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
                double rj = (n_psi-1.0)*log(rj1)/log(rj2);

                int j = int(rj);

                if (j < 0)
                    j = 0;
                else if (j >= n_psi-1)
                    j = n_psi - 2;

                double frac = rj - j;

                double bulge_dens_psi = *(Dens_Psi_Bulge_CUDA+j) + 
                                       frac*(*(Dens_Psi_Bulge_CUDA+j+1)-*(Dens_Psi_Bulge_CUDA+j));
                totdens += bulge_dens_psi;
            }
            
            if (j==0 || j==n_theta)
            else if (j%2==0)
            {
                totdens *= 2;
            }
            else if (j%2==1)
            {
                totdens *= 4;
            }
        }
        
        totdens *= d_cos_theta*fourpi*0.3333333333333;
        return totdens;
    }
};

__device__ double PotCUDA(double s, double z, thrust::device_ptr<double> A_Pot)
{ 
    double pot = 0, r=sqrt(s*s+z*z);
    
    if (r==0)
    {
        return a00*oneoversqrt4pi;
    }
    
    int ihi = ceil(r/dr);
    
    if (r<dr)
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
        //cout << "Pot finds out of range indices. Possible int overflow. " <<
        //        "Exiting..." << endl;
        //cout <<"ihi = "<<ihi<<" "<<r<<" "<<dr<<" "<<r/dr<<" "<<ceil(r/dr)<<" "
        //     <<s<<" "<<z<<endl;
        exit(1);
    }
    else if (ihi > nr-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<s<<" "<<z<<" "<< endl;
        ihi = nr-1;
    }
    
    double r1 = *(Radius+ihi-1);//Radius[ihi-1];
    double r2 = *(Radius+ihi);//Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    double cos_theta = z/r;
    
    pot += LegendreCUDA(l, cos_theta)*Plcon[l]*
           (t*(*(A_Pot+ihi)) + tm1*A_Pot(*(A_Pot+ihi-1)));
        
    return pot;
}

__device__ double LegendreCUDA(const int l, const double x)
{
    if (l<0 || x<-1 || x>1)
    {
        cout << "Error! out of range legendre arguments" << endl;
        exit(1);
    }
    else if (l==0)
    {
        return 1;
    }
    else if (l==1)
    {
        return x;
    }
    else if (l==2)
    {
        return 0.5*(3*x*x-1);
    }
    else if (x = 1.0)
    {
        return 1;
    }
    else if (x==-1)
    {
        if (l%2==1)
        {
            return -1;
        }
        else if (l%2==0)
        {
            return 1;
        }
    }
    else
    {
        double p_ellm2 = 1;
        double p_ellm1 = x;
        double p_ell = p_ellm1;
        
        //double p_ellm2 = EPSILON;
        
        for (int ell=2; ell<l+1; ++ell)
        {
            p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2)/ell;
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
        }
        
        return p_ell;
    }
}

double GetADensCUDA(void)
{
    cudaMemcpyToSymbol(halo_flag_d, &halo_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_flag_d, &disk_flag, sizeof(int));
    cudaMemcpyToSymbol(bulge_flag_d, &bulge_flag, sizeof(int));
    cudaMemcpyToSymbol(n_psi_d, &n_psi, sizeof(int));
    cudaMemcpyToSymbol(nr_d, &nr, sizeof(int));
    cudaMemcpyToSymbol(lmax_d, &lmax, sizeof(int));
    cudaMemcpyToSymbol(l_max_d, &l_max, sizeof(int));
    cudaMemcpyToSymbol(gasdisk_flag_d, &gasdisk_flag, sizeof(int));
    cudaMemcpyToSymbol(smbh_flag_d, &smbh_flag, sizeof(int));
    cudaMemcpyToSymbol(nbody_flag_d, &nbody_flag, sizeof(int));
    cudaMemcpyToSymbol(do_file_io_d, &do_file_io, sizeof(int));
    cudaMemcpyToSymbol(chisq_flag_d, &chisq_flag, sizeof(int));
    cudaMemcpyToSymbol(sersic_flag_d, &sersic_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_d, &disk, sizeof(int));
    cudaMemcpyToSymbol(gas_disk_d, &gas_disk, sizeof(int));
    cudaMemcpyToSymbol(gasdisk_params_d, &gasdisk_params, sizeof(int));
    cudaMemcpyToSymbol(disk_params_d, &disk_params, sizeof(int));
    cudaMemcpyToSymbol(nondisk_params_d, &nondisk_params, sizeof(int));
    cudaMemcpyToSymbol(astro_params_d, &astro_params, sizeof(int));
    cudaMemcpyToSymbol(error_params_d, &error_params, sizeof(int));
    cudaMemcpyToSymbol(dr_d, &dr, sizeof(int));
    cudaMemcpyToSymbol(r_max_d, &r_max, sizeof(int));
    cudaMemcpyToSymbol(n_int_d, &n_int, sizeof(int));
    cudaMemcpyToSymbol(nr_spline_d, &nr_spline, sizeof(int));
    cudaMemcpyToSymbol(n_iter_d, &n_iter, sizeof(int));
    cudaMemcpyToSymbol(max_iter_d, &max_iter, sizeof(int));
    cudaMemcpyToSymbol(iter_d, &iter, sizeof(int));
    cudaMemcpyToSymbol(n_simpson_d, &n_simpson, sizeof(int));
    cudaMemcpyToSymbol(dr_d, &dr, sizeof(double));
    cudaMemcpyToSymbol(r_max_d, &r_max, sizeof(double));
    cudaMemcpyToSymbol(potcor_d, &potcor, sizeof(double));
    cudaMemcpyToSymbol(psi_crit_d, &psi_crit, sizeof(double));
    cudaMemcpyToSymbol(psi_d_d, &psi_d, sizeof(double));
    cudaMemcpyToSymbol(psi_0_d, &psi_0, sizeof(double));
    cudaMemcpyToSymbol(tolerance_factor_d, &tolerance_factor, sizeof(double));
    cudaMemcpyToSymbol(fraction_d, &fraction, sizeof(double));
    cudaMemcpyToSymbol(psi_upper_initial_d, &psi_upper_initial, sizeof(double));
    cudaMemcpyToSymbol(log_rmax_d, &log_rmax, sizeof(double));
    cudaMemcpyToSymbol(log_dr_d, &log_dr, sizeof(double));
    cudaMemcpyToSymbol(delta_logr_d, &delta_logr, sizeof(double));
    cudaMemcpyToSymbol(total_mass_d, &total_mass, sizeof(double));
    cudaMemcpyToSymbol(halo_mass_d, &halo_mass, sizeof(double));
    cudaMemcpyToSymbol(total_disk_mass_d, &total_disk_mass, sizeof(double));
    cudaMemcpyToSymbol(bulge_mass_d, &bulge_mass, sizeof(double));
    cudaMemcpyToSymbol(r_edge_d, &r_edge, sizeof(double));
    cudaMemcpyToSymbol(halo_edge_d, &halo_edge, sizeof(double));
    cudaMemcpyToSymbol(disk_edge_d, &disk_edge, sizeof(double));
    cudaMemcpyToSymbol(bulge_edge_d, &bulge_edge, sizeof(double));
    cudaMemcpyToSymbol(Q_total_d, &Q_total, sizeof(double));
    cudaMemcpyToSymbol(sin_inclination_d, &sin_inclination, sizeof(double));
    cudaMemcpyToSymbol(cos_inclination_d, &cos_inclination, sizeof(double));
    cudaMemcpyToSymbol(Plcon_d, Plcon, sizeof(double));
    
    thrust::device_ptr<float> A_Pot_CUDA = thrust::device_malloc<float>(nr);
    thrust::device_ptr<float> A_Dens_CUDA = thrust::device_malloc<float>(nr);
    thrust::device_ptr<float> Disk_Density = thrust::device_malloc<float>(nr);
    thrust::device_ptr<float> Halo_Density = thrust::device_malloc<float>(nr);
    thrust::device_ptr<float> Bulge_Density = thrust::device_malloc<float>(nr);

    thrust::copy(A_Dens[l/2].start(), A_Dens[l/2].end(), A_Dens_CUDA);
    thrust::copy(A_Pot[l/2].start(), A_Pot[l/2].end(), A_Pot_CUDA);

    thrust::transform(radius_dev_vtr,Disk_Density.start(),GetDiskDensity(G));
    thrust::transform(radius_dev_vtr,Halo_Density.start(),GetHaloDensity(G));
    thrust::transform(radius_dev_vtr,Bulge_Density.start(),GetBulgeDensity(G));
    //potential should be a device pointer
    thrust::transform(Radius.begin(), Radius.end(), A_Dens_CUDA.begin(), GetDensity());
    
    thrust::copy(A_Dens_CUDA.start(), A_Dens_CUDA.end(), A_Dens[l/2]);
}
