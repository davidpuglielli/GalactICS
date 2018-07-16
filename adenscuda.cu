#include "galaxy_cuda.h"

void CopyGlobalsToDevice1(void)
{
    cudaMemcpyToSymbol(c_halo_d, &G.c_halo, sizeof(G.c_halo));
    cudaMemcpyToSymbol(v_halo_d, &G.v_halo, sizeof(G.v_halo));
    cudaMemcpyToSymbol(a_halo_d, &G.a_halo, sizeof(G.a_halo));
    cudaMemcpyToSymbol(drtrunc_halo_d, &G.drtrunc_halo, sizeof(G.drtrunc_halo));
    cudaMemcpyToSymbol(cusp_d, &G.cusp, sizeof(G.cusp));
    cudaMemcpyToSymbol(halo_stream_d, &G.halo_stream, sizeof(G.halo_stream));
    cudaMemcpyToSymbol(n_sersic_d, &G.n_sersic, sizeof(G.n_sersic));
    cudaMemcpyToSymbol(ppp_d, &G.ppp, sizeof(G.ppp));
    cudaMemcpyToSymbol(b_n_d, &G.b_n, sizeof(G.b_n));
    cudaMemcpyToSymbol(v_bulge_d, &G.v_bulge, sizeof(G.v_bulge));
    cudaMemcpyToSymbol(a_bulge_d, &G.a_bulge, sizeof(G.a_bulge));
    cudaMemcpyToSymbol(bulge_stream_d, &G.bulge_stream, sizeof(G.bulge_stream));
    cudaMemcpyToSymbol(rho_0_d, &G.rho_0, sizeof(G.rho_0));
    cudaMemcpyToSymbol(bh_mass_d, &G.bh_mass, sizeof(G.bh_mass));
    
    cudaMemcpyToSymbol(M_Disk, &G.M_Disk[0], sizeof(G.M_Disk));
    cudaMemcpyToSymbol(R_Disk, &G.R_Disk[0], sizeof(G.R_Disk));
    cudaMemcpyToSymbol(Z_Disk, &G.Z_Disk[0], sizeof(G.Z_Disk));
    cudaMemcpyToSymbol(Out_Disk, &G.Out_Disk[0], sizeof(G.Out_Disk));
    cudaMemcpyToSymbol(Dr_Trunc, &G.Dr_Trunc[0], sizeof(G.Dr_Trunc));
    if (!G.R_Kormendy.empty() && !G.Alpha.empty())
    {        
        cudaMemcpyToSymbol(R_Kormendy, &G.R_Kormendy[0], sizeof(G.R_Kormendy));
        cudaMemcpyToSymbol(Alpha, &G.Alpha[0], sizeof(G.Alpha));
    }
    cudaMemcpyToSymbol(Sigma_0, &G.Sigma_0[0], sizeof(G.Sigma_0));
    cudaMemcpyToSymbol(R_Sigma, &G.R_Sigma[0], sizeof(G.R_Sigma));
    cudaMemcpyToSymbol(Disk_Const_d, &Disk_Const[0], sizeof(Disk_Const));
    cudaMemcpyToSymbol(Rho_Disk_Const_d, &Rho_Disk_Const[0], sizeof(Rho_Disk_Const));
    
    if (gasdisk_flag)
    {
        cudaMemcpyToSymbol(Z_GasDisk, &G.Z_GasDisk[0], sizeof(G.Z_GasDisk));
        cudaMemcpyToSymbol(R_GasDisk, &G.R_GasDisk[0], sizeof(G.R_GasDisk));
        cudaMemcpyToSymbol(Z_GasDisk, &G.Z_GasDisk[0], sizeof(G.Z_GasDisk));
        cudaMemcpyToSymbol(Out_GasDisk, &G.Out_GasDisk[0], sizeof(G.Out_GasDisk));
        cudaMemcpyToSymbol(Dr_Trunc_Gas, &G.Dr_Trunc_Gas[0], sizeof(G.Dr_Trunc_Gas));
        if (!G.R_Kormendy_Gas.empty() && !G.Alpha_Gas.empty())
        {        
            cudaMemcpyToSymbol(R_Kormendy_Gas, &G.R_Kormendy_Gas[0], sizeof(G.R_Kormendy_Gas));
            cudaMemcpyToSymbol(Alpha_Gas, &G.Alpha_Gas[0], sizeof(G.Alpha_Gas));
        }
        cudaMemcpyToSymbol(Sigma_0_Gas, &G.Sigma_0_Gas[0], sizeof(G.Sigma_0_Gas));
        cudaMemcpyToSymbol(R_Sigma_Gas, &G.R_Sigma_Gas[0], sizeof(G.R_Sigma_Gas));
        cudaMemcpyToSymbol(GasDisk_Const_d, &GasDisk_Const[0], sizeof(GasDisk_Const));
        cudaMemcpyToSymbol(Rho_GasDisk_Const_d, &Rho_GasDisk_Const[0], sizeof(Rho_GasDisk_Const));
    }
    
    cudaMemcpyToSymbol(halo_flag_d, &halo_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_flag_d, &disk_flag, sizeof(int));
    cudaMemcpyToSymbol(bulge_flag_d, &bulge_flag, sizeof(int));
    cudaMemcpyToSymbol(n_psi_d, &n_psi, sizeof(int));
    cudaMemcpyToSymbol(nr_d, &nr, sizeof(int));
    cudaMemcpyToSymbol(smbh_flag_d, &smbh_flag, sizeof(int));
    cudaMemcpyToSymbol(nbody_flag_d, &nbody_flag, sizeof(int));
    cudaMemcpyToSymbol(do_file_io_d, &do_file_io, sizeof(int));
    cudaMemcpyToSymbol(chisq_flag_d, &chisq_flag, sizeof(int));
    cudaMemcpyToSymbol(sersic_flag_d, &sersic_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_d, &disk, sizeof(int));
    
    cudaMemcpyToSymbol(l_max_d, &l_max, sizeof(int));
    cudaMemcpyToSymbol(n_int_d, &n_int, sizeof(int));
    cudaMemcpyToSymbol(nr_spline_d, &nr_spline, sizeof(int));
    cudaMemcpyToSymbol(n_iter_d, &n_iter, sizeof(int));
    cudaMemcpyToSymbol(max_iter_d, &max_iter, sizeof(int));
    cudaMemcpyToSymbol(iter_d, &iter, sizeof(int));
    cudaMemcpyToSymbol(n_simpson_d, &n_simpson, sizeof(int));
    
    cudaMemcpyToSymbol(dr_d, &dr, sizeof(double));
    cudaMemcpyToSymbol(r_max_d, &r_max, sizeof(double));
    cudaMemcpyToSymbol(tolerance_factor_d, &tolerance_factor, sizeof(double));
    cudaMemcpyToSymbol(fraction_d, &fraction, sizeof(double));
    cudaMemcpyToSymbol(log_rmax_d, &log_rmax, sizeof(double));
    cudaMemcpyToSymbol(log_dr_d, &log_dr, sizeof(double));
    cudaMemcpyToSymbol(delta_logr_d, &delta_logr, sizeof(double));
    
    cudaMemcpyToSymbol(potcor_d, &potcor, sizeof(double));
    cudaMemcpyToSymbol(psi_crit_d, &psi_crit, sizeof(double));
    cudaMemcpyToSymbol(psi_d_d, &psi_d, sizeof(double));
    cudaMemcpyToSymbol(psi_0_d, &psi_0, sizeof(double));
    cudaMemcpyToSymbol(log_rj2_d, &log_rj2, sizeof(double));
    cudaMemcpyToSymbol(psi_upper_initial_d, &psi_upper_initial, sizeof(double));
    cudaMemcpyToSymbol(Plcon_d, Plcon, 40*sizeof(double));
}
    
struct TotalDensityFunctor
{
    const thrust::device_ptr<double> Radius_CUDA;
    const thrust::device_ptr<double> A_Pot_CUDA;
    const thrust::device_ptr<double> Dens_Psi_Halo_CUDA;
    const thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
    
    TotalDensityFunctor(thrust::device_ptr<double> _Radius,
                        thrust::device_ptr<double> _A_Pot, 
                        thrust::device_ptr<double> _Dens_Psi_Halo,
                        thrust::device_ptr<double> _Dens_Psi_Bulge):
                        Radius_CUDA(_Radius), A_Pot_CUDA(_A_Pot), 
                        Dens_Psi_Halo_CUDA(_Dens_Psi_Halo), 
                        Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge){}
    //saxpy_functor(galstruct _a, double _b):a(_a), b(_b){}
    TotalDensityFunctor(){}
    
    __host__ __device__ 
    float operator()(const float& rad) const
    {
        //printf("Now %f\n", rad);
        
        int n_theta = 20;
        double d_cos_theta = 1.0/n_theta;
        
        double dens_harmonic = 0;
        double disk_dens, halo_dens_psi, bulge_dens_psi, appdiskdens;
            
        for (int j = 0; j < n_theta+1; ++j)
        {
            //printf("Now %f\n", rad);
            //printf("index %d %f\n", j, rad);
            
            double cos_theta = j*d_cos_theta;
            double z = rad*cos_theta;
            double ss = rad*sqrt(1-cos_theta*cos_theta);
            
            double totdens = 0;

            double psi = PotCUDA(ss, z, Radius_CUDA, A_Pot_CUDA);
            
            disk_dens = 0; 
            halo_dens_psi = 0; 
            bulge_dens_psi = 0;
            appdiskdens = 0;
            
            //printf("psi  %d %f %f %f    ", j, rad, z, psi);

            for (int i = 0; i < disk_d; ++i)
            {
                if(fabs(z/(*(Z_Disk+i)))>30)
                    continue;

                double trunc_fac = GetTruncCUDA(ss, i);

                if (trunc_fac==0)
                    continue;

                double con;

                if (z==0)
                {
                    con = 1;
                }
                else if (!halo_flag_d && !bulge_flag_d)
                {
                    con = exp(-fabs(z/(*(Z_Disk+i))));
                    con = 2.0*con/(1.0+con*con);
                    con *= con;
                }
                else
                {
                    //KK: make disk vertically isothermal; normalize density so that
                    //at height 3*zdisk the vertical density drop is that of a sech^2.
                    double psizh = PotCUDA(ss, 3.0*(*(Z_Disk+i)), 
                                           Radius_CUDA, A_Pot_CUDA);
                    double psi00 = PotCUDA(ss, 0, Radius_CUDA, A_Pot_CUDA);
                    double dpsizh = psi00 - psizh;
                    double dpsi = psi00 - psi;
                    double coeff = 0;

                    if (ss==0 && z==0)
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

                DiskProfileCUDA(ss, z, i, Sigma_Profile, Rho_Profile);
                disk_dens += (*(Rho_Disk_Const_d+i))*Sigma_Profile[0]*con*trunc_fac;
                //printf("dens  %d %f %f %f %f  %f    %e %e %e %e\n", j, ss, rad, 
                //        sqrt(1-cos_theta*cos_theta), z, psi, 
                //        *(Rho_Disk_Const_d+i),Sigma_Profile[0],con,trunc_fac);
                //totdens += disk_dens;
                //printf("  %f ", totdens);
            }
            
            totdens += disk_dens;

            //Do the halo            
            if (psi < psi_crit_d)
            {
                totdens+=0;
                halo_dens_psi = 0;
            }                    
            else if (psi >= psi_0_d) 
            {
                totdens += *(Dens_Psi_Halo_CUDA);
                halo_dens_psi = *(Dens_Psi_Halo_CUDA);
            }
            else
            {
                double rj1 = (psi_0_d - psi)/(psi_0_d - psi_d_d);
                double rj = (n_psi_d-1.0)*log(rj1)*log_rj2_d;

                int j = int(rj);

                if (j < 0)
                    j = 0;
                else if (j >= n_psi_d-1)
                    j = n_psi_d - 2;

                double frac = rj - j;

                halo_dens_psi = *(Dens_Psi_Halo_CUDA+j) + 
                                frac*(*(Dens_Psi_Halo_CUDA+j+1)-(*(Dens_Psi_Halo_CUDA+j)));
                totdens += halo_dens_psi;
            }

            //Do the bulge
            if (psi < psi_crit_d)
            {
                totdens+=0;
                bulge_dens_psi = 0;
            }                    
            else if (psi >= psi_0_d) 
            {
                totdens += *(Dens_Psi_Bulge_CUDA);
                bulge_dens_psi = *(Dens_Psi_Bulge_CUDA);
            }
            else
            {
                double rj1 = (psi_0_d - psi)/(psi_0_d - psi_d_d);
                double rj = (n_psi_d-1.0)*log(rj1)*log_rj2_d;

                int j = int(rj);

                if (j < 0)
                    j = 0;
                else if (j >= n_psi_d-1)
                    j = n_psi_d - 2;

                double frac = rj - j;

                bulge_dens_psi = *(Dens_Psi_Bulge_CUDA+j) + 
                                 frac*(*(Dens_Psi_Bulge_CUDA+j+1)-*(Dens_Psi_Bulge_CUDA+j));
                totdens += bulge_dens_psi;
            }
            
            //Do the high-frequency disk components
            for (int i = 0; i < disk_d; ++i)
            {
                double trunc_fac, trunc_facprime;

                GetTruncPrimeCUDA(rad, i, trunc_fac, trunc_facprime);

                double Sigma_Profile[3], Rho_Profile[3];

                DiskProfile2PrimeCUDA(rad, z, i, Sigma_Profile, Rho_Profile);

                double f = Sigma_Profile[0]*(*(Z_Disk+i))*(*(Z_Disk+i))*trunc_fac;
                double f1r = (*(Z_Disk+i))*(*(Z_Disk+i))*
                             (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/rad;
                double f2 = (*(Z_Disk+i))*(*(Z_Disk+i))*(Sigma_Profile[2]*trunc_fac +
                            2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		                    (rad-(*(Out_Disk+i)))/(*(Dr_Trunc+i))/(*(Dr_Trunc+i)));

                if (rad==0)
                {
                    f1r = f2 = 0;
                }        

                double zz = fabs(z/(*(Z_Disk+i)));
                double ezz = exp(-zz), e2zz = exp(-2*zz);
                double tlncosh = zz+log(0.5*(1.0 + e2zz));
                double tztanh = zz*(1-e2zz)/(1+e2zz);
                double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
                double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
                               f*tsech2/(*(Z_Disk+i))/(*(Z_Disk+i));
                appdiskdens += (*(Rho_Disk_Const_d+i))*total;
            }
            
            totdens -= appdiskdens;
            //printf("psi  %d %f %f %f     %e %e %e %e\n", j, rad, z, psi, 
            //        disk_dens, halo_dens_psi, bulge_dens_psi, appdiskdens);
            totdens *= LegendreCUDA(l_d, cos_theta)*Plcon_d[l_d];
            
            if (j==0 || j==n_theta)
                ;
            else if (j%2==0)
            {
                totdens *= 2;
            }
            else if (j%2==1)
            {
                totdens *= 4;
            }
            
            dens_harmonic += totdens;
            //printf("In loop %e %e %e   %e %e %e %e     %e %e %e\n", rad, ss, z, 
            //        disk_dens, halo_dens_psi, bulge_dens_psi, appdiskdens, 
            //        totdens, dens_harmonic, LegendreCUDA(l_d, cos_theta));
        }
        
        //printf("dens_harmonic %f %f %f %f %f\n",dens_harmonic, rad, psi_0_d, psi_crit_d, psi_d_d);
        //printf("Now %f\n", rad);//*A_Pot_0_CUDA[3]);
        
        dens_harmonic *= d_cos_theta*fourpi*0.3333333333333;
        
        //printf("In CUDA %e %e %e %e %e %e\n", rad, dens_harmonic, disk_dens, 
        //       halo_dens_psi, bulge_dens_psi, appdiskdens);
        return dens_harmonic;
    }
};

void GetADensCUDA(const int l)
{
    cudaMemcpyToSymbol(l_d, &l, sizeof(int));
    
    cudaMemcpyToSymbol(lmax_d, &lmax, sizeof(int));
    cudaMemcpyToSymbol(Q_total_d, &Q_total, sizeof(double));
    cudaMemcpyToSymbol(a00_d, &A_Pot[0][0], sizeof(double));
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2);
    thrust::device_ptr<double> A_Dens_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Halo_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA = thrust::device_malloc<double>(nr);

    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    thrust::copy(A_Dens[l/2].begin(), A_Dens[l/2].end(), A_Dens_CUDA);
    
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[l/2].begin(), A_Pot[l/2].end(), A_Pot_CUDA+offset);
    }
    
    thrust::copy(Dens_Psi_Halo.begin(), Dens_Psi_Halo.end(), Dens_Psi_Halo_CUDA);
    thrust::copy(Dens_Psi_Bulge.begin(), Dens_Psi_Bulge.end(), Dens_Psi_Bulge_CUDA);

    thrust::transform(Radius_CUDA, Radius_CUDA+nr, A_Dens_CUDA, 
                      TotalDensityFunctor(Radius_CUDA, A_Pot_CUDA,
                                          Dens_Psi_Halo_CUDA, Dens_Psi_Bulge_CUDA));
    //thrust::transform(Radius_CUDA, Radius_CUDA+nr, Dens_Psi_Halo_CUDA, A_Dens_CUDA, 
    //                  TotalDensityFunctor());
    //cout << nr_d << " " << Radius_CUDA[100] << " " << A_Pot_CUDA[100] << " " << A_Dens_CUDA[100] << endl;
    
    thrust::copy(A_Dens_CUDA, A_Dens_CUDA+nr, A_Dens[l/2].begin());
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(A_Dens_CUDA);
    thrust::device_free(Dens_Psi_Halo_CUDA);
    thrust::device_free(Dens_Psi_Bulge_CUDA);
}
__host__ __device__ 
double PotCUDA(double s, double z, 
               thrust::device_ptr<double> Radius_CUDA,
               thrust::device_ptr<double> A_Pot_CUDA)
{ 
    double pot = 0, r=sqrt(s*s+z*z);
    
    if (r==0)
    {
        return a00_d*oneoversqrt4pi;
    }
    
    int ihi = ceil(r/dr_d);
    
    if (r<dr_d)
    {
        ihi = 1;
    }
    else if (r > 10*r_max_d)
    {
        ihi = nr_d-1;
    }
    else if (ihi < 1)
    {
    }
    else if (ihi > nr_d-1)
    {
        ihi = nr_d-1;
    }
    
    double r1 = *(Radius_CUDA+ihi-1);
    double r2 = *(Radius_CUDA+ihi);
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    double cos_theta = z/r;
    
    for (int l = l_max_d; l > -1; l-=2)
    {
        pot += LegendreCUDA(l, cos_theta)*Plcon_d[l]*
               (t*(*(A_Pot_CUDA+l/2*nr_d+ihi)) + tm1*(*(A_Pot_CUDA+l/2*nr_d+ihi-1)));
    }
    
//     pot += LegendreCUDA(0, cos_theta)*Plcon_d[0]*
//            (t*(*(A_Pot_0_CUDA+ihi)) + tm1*(*(A_Pot_0_CUDA+ihi-1)));
//         
//     pot += LegendreCUDA(2, cos_theta)*Plcon_d[2]*
//            (t*(*(A_Pot_2_CUDA+ihi)) + tm1*(*(A_Pot_2_CUDA+ihi-1)));
        
    if (disk_flag_d)
    {
        for (int i = 0; i < disk_d; ++i)
        {
            double Sigma_Profile[2], Rho_Profile[2];

            DiskProfileCUDA(r, z, i, Sigma_Profile, Rho_Profile);

            double trunc_fac = GetTruncCUDA(r, i);
            double zz = fabs(z/(*(Z_Disk+i)));
            pot += -fourpi*(*(Rho_Disk_Const_d+i))*Sigma_Profile[0]*(*(Z_Disk+i))*
                   (*(Z_Disk+i))*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
        }
    }
    
    //printf("psir %f %f      %f %f %f      %f %f %f\n", 
    //        s, z, 1.0*(*(A_Pot_0_CUDA+ihi)), 1.0*(*(A_Pot_0_CUDA+ihi-1)), t,
    //        1.0*(*(A_Pot_2_CUDA+ihi)), 1.0*(*(A_Pot_2_CUDA+ihi-1)), tm1);
    
    return pot;
}

__host__ __device__ 
double LegendreCUDA(const int l_d, const double x)
{
    if (l_d<0 || x<-1 || x>1)
    {
        //cout << "Error! out of range legendre arguments" << endl;
        //exit(1);
    }
    else if (l_d==0)
    {
        return 1;
    }
    else if (l_d==1)
    {
        return x;
    }
    else if (l_d==2)
    {
        return 0.5*(3*x*x-1);
    }
    else if (x==1.0)
    {
        return 1;
    }
    else if (x==-1)
    {
        if (l_d%2==1)
        {
            return -1;
        }
        else if (l_d%2==0)
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
        
        for (int ell=2; ell<l_d+1; ++ell)
        {
            p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2)/ell;
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
        }
        
        return p_ell;
    }
}

__host__ __device__ 
void DiskProfileCUDA(const double radius, const double z, const int i, 
                     double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0;//, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_Disk+i)));
        //ddensity = -density / G.(*(R_Disk+i));
        //dddensity = -ddensity / G.(*(R_Disk+i));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
                      cosh(fabs(z/(*(Z_Disk+i)))));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_Disk+i)));
        //ddensity = -density / G.(*(R_Disk+i));
        //dddensity = -ddensity / G.(*(R_Disk+i));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
                      cosh(fabs(z/(*(Z_Disk+i)))));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    //Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}

__host__ __device__ 
void DiskProfile2PrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_Disk+i)));
        ddensity = -density / (*(R_Disk+i));
        dddensity = -ddensity / (*(R_Disk+i));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
                      cosh(fabs(z/(*(Z_Disk+i)))));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_Disk+i)));
        ddensity = -density / (*(R_Disk+i));
        dddensity = -ddensity / (*(R_Disk+i));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
                      cosh(fabs(z/(*(Z_Disk+i)))));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    Sigma_Profile[1] = ddensity;
    Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}

__host__ __device__ 
double GetTruncCUDA(const double &rad, const int &i)
{
    //double truncfac;
    double erfarg = (rad - (*(Out_Disk+i))) * oneoversqrt2 / (*(Dr_Trunc+i));
    
    if (erfarg > 4) return 0;
    else if (erfarg < -4) return 1;
    else return 0.5 * erfc(erfarg);
}

__host__ __device__ 
void GetTruncPrimeCUDA(const double &r, const int &i, 
                       double &truncfac, double &truncfacprime)
{
    double erfarg = (r - (*(Out_Disk+i)))*oneoversqrt2/(*(Dr_Trunc+i));
    
    if (erfarg > 4) 
    {
        truncfac = 0;
        truncfacprime = 0;
    }
    else if (erfarg < -4) 
    {
        truncfac = 1;
        truncfacprime = 0;
    }
    else
    {
        truncfac = 0.5 * erfc(erfarg);
        truncfacprime = -exp(-erfarg*erfarg)/(*(Dr_Trunc+i))*oneoversqrt2pi;
    }
}

