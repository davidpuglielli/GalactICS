struct DiskDensityEstimateFunctor
{
    //TotalDensityFunctor(){}
    
    __host__ __device__ 
    float operator()(const float& rad) const
    {
        int n_theta = 100;
        double cos_theta_max = min(1.0, 10*Z_Disk[0]/rad);
        double d_cos_theta = cos_theta_max/n_theta;
        
        double density = 0;
            
        for (int j = 0; j < n_theta+1; ++j)
        {
            double cos_theta = j*d_cos_theta;
            double z = rad*cos_theta;
            double ss = rad*sqrt(1.0-cos_theta*cos_theta);
            double dens = 0;
            
            for (int i = 0; i < disk_d; ++i)
            {
                double trunc_fac = GetTruncCUDA(rad, i);
	            double Sigma_Profile[2], Rho_Profile[2];

                DiskProfileCUDA(rad, z, i, Sigma_Profile, Rho_Profile);
                dens += (*(Rho_Disk_Const_d+i))*Rho_Profile[0]*trunc_fac;
            }
            
            if (j==0 || j==n_theta)
                ;
            else if (j%2==0)
            {
                dens *= 2;
            }
            else if (j%2==1)
            {
                dens *= 4;
            }
            
            density += dens;
        }
        
        density *= d_cos_theta*0.3333333333333;
        
        return density;
    }
};

struct TotalDensityFunctor
{
    const thrust::device_ptr<double> Radius_CUDA;
    const thrust::device_ptr<double> A_Pot_CUDA;
    const thrust::device_ptr<double> Dens_Psi_Halo_CUDA;
    const thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
    const thrust::device_ptr<double> GasDensity_Const_CUDA;
    const thrust::device_ptr<double> Polytrope_Const_CUDA;
    
    TotalDensityFunctor(thrust::device_ptr<double> _Radius,
                        thrust::device_ptr<double> _A_Pot, 
                        thrust::device_ptr<double> _Dens_Psi_Halo,
                        thrust::device_ptr<double> _Dens_Psi_Bulge,
                        thrust::device_ptr<double> _GasDensity_Const_CUDA,
                        thrust::device_ptr<double> _Polytrope_Const_CUDA):
                        Radius_CUDA(_Radius), A_Pot_CUDA(_A_Pot), 
                        Dens_Psi_Halo_CUDA(_Dens_Psi_Halo), 
                        Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge),
                        GasDensity_Const_CUDA(_GasDensity_Const_CUDA),
                        Polytrope_Const_CUDA(_Polytrope_Const_CUDA){}
    //saxpy_functor(galstruct _a, double _b):a(_a), b(_b){}
    TotalDensityFunctor(){}
    
    __host__ __device__ //__noinline__ 
    float operator()(const float& rad) const
    {
        //printf("Now %f\n", rad);
        
        int n_theta = 100;
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
            //printf("Now %f\n", rad);
            double psi = PotCUDA(ss, z, Radius_CUDA, A_Pot_CUDA);
            //printf("Now1 %f\n", rad);            
            disk_dens = 0; 
            halo_dens_psi = 0; 
            bulge_dens_psi = 0;
            appdiskdens = 0;
            
            //printf("psi  %d %f %f %f    ", j, rad, z, psi);

            //Do the disk(s)
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
                
                //con *= (1-exp(-ss*ss/5));
//     double zfac = (*(Z_Disk+i))*(1-0.95*exp(-ss*ss/25));//+0.00001;
//     con = 1.0/cosh(fabs(z/zfac))/cosh(fabs(z/zfac));

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

        //printf("Now %f\n", rad);
        
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

        //printf("Now %f\n", rad);
        
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
            
        //printf("Now %f\n", rad);
        
            //Do the high-frequency disk components
            appdiskdens = AppDiskDensCUDA(ss, z);
//             for (int i = 0; i < disk_d; ++i)
//             {
//                 double trunc_fac, trunc_facprime;
// 
//                 GetTruncPrimeCUDA(rad, i, trunc_fac, trunc_facprime);
// 
//                 double Sigma_Profile[3], Rho_Profile[3];
// 
//                 DiskProfile2PrimeCUDA(rad, z, i, Sigma_Profile, Rho_Profile);
// 
//                 double f = Sigma_Profile[0]*(*(Z_Disk+i))*(*(Z_Disk+i))*trunc_fac;
//                 double f1r = (*(Z_Disk+i))*(*(Z_Disk+i))*
//                              (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/rad;
//                 double f2 = (*(Z_Disk+i))*(*(Z_Disk+i))*(Sigma_Profile[2]*trunc_fac +
//                             2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
// 		                    (rad-(*(Out_Disk+i)))/(*(Dr_Trunc+i))/(*(Dr_Trunc+i)));
// 
//                 if (rad==0)
//                 {
//                     f1r = f2 = 0;
//                 }        
// 
//                 double zz = fabs(z/(*(Z_Disk+i)));
//                 double ezz = exp(-zz), e2zz = exp(-2*zz);
//                 double tlncosh = zz+log(0.5*(1.0 + e2zz));
//                 double tztanh = zz*(1-e2zz)/(1+e2zz);
//                 double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
//                 double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
//                                f*tsech2/(*(Z_Disk+i))/(*(Z_Disk+i));
//                 appdiskdens += (*(Rho_Disk_Const_d+i))*total;
//             }
            
            totdens -= appdiskdens;
            //printf("psi  %d %f %f %f     %e %e %e %e\n", j, rad, z, psi, 
            //        disk_dens, halo_dens_psi, bulge_dens_psi, appdiskdens);
            
            //Do the gas disk(s)
            double gas_dens = 0;
            for (int i = 0; i < gas_disk_d; ++i)
            {
                gas_dens += GasDiskDensfICUDA(ss, z, i, Radius_CUDA,
                                              A_Pot_CUDA, Polytrope_Const_CUDA,
                                              GasDensity_Const_CUDA);
            }
            
            totdens += gas_dens;
            
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
        
       // printf("Now %f\n", rad);
        
        //printf("In CUDA %e %e %e %e %e %e\n", rad, dens_harmonic, disk_dens, 
        //       halo_dens_psi, bulge_dens_psi, appdiskdens);
        return dens_harmonic;
    }
};
