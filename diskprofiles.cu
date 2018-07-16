__host__ __device__     
void DiskProfileCUDA(const double radius, const double z, const int i, 
                     double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0;//, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_Disk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_Disk+i)));
        //ddensity = -density / G.(*(R_Disk+i));
        //dddensity = -ddensity / G.(*(R_Disk+i));
        //rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
        //              cosh(fabs(z/(*(Z_Disk+i)))));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
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
void DiskProfilePrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0;//, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_Disk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_Disk+i)));
        ddensity = -density / (*(R_Disk+i));
        //dddensity = -ddensity / (*(R_Disk+i));
        //rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
        //              cosh(fabs(z/(*(Z_Disk+i)))));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_Disk+i)));
        ddensity = -density / (*(R_Disk+i));
        //dddensity = -ddensity / (*(R_Disk+i));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
                      cosh(fabs(z/(*(Z_Disk+i)))));
    }
    
    //Inner Truncated exponential
    
    //Thick disk
    
    //etc.
    
    Sigma_Profile[0] = density;
    Sigma_Profile[1] = ddensity;
    //Sigma_Profile[2] = dddensity;
    
    Rho_Profile[0] = rho_density;
}
  
__host__ __device__     
void DiskProfile2PrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_Disk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_Disk+i)));
        ddensity = -density / (*(R_Disk+i));
        dddensity = -ddensity / (*(R_Disk+i));
        //rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/(*(Z_Disk+i))))*
        //              cosh(fabs(z/(*(Z_Disk+i)))));
        rho_density = exp(-radius/(*(R_Disk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
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
void AppDiskForceCUDA(double &r, double &z, double &fsad, double &fzad)
{
    double s = sqrt(r*r + z*z);
    fsad = fzad = 0;
    //return;
        
    for (int i = 0; i < disk_d; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrimeCUDA(s, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[2], Rho_Profile[2];

        DiskProfilePrimeCUDA(s, z, i, Sigma_Profile, Rho_Profile);
    
        double zz = fabs(z/(*(Z_Disk+i)));
        double e2zz = exp(-2*zz);
        double tlncoshz = zz+log(0.5*(1.0 + e2zz));
        
		double f = Sigma_Profile[0]*trunc_fac;
		double f_prime = Sigma_Profile[0]*trunc_facprime + 
                         Sigma_Profile[1]*trunc_fac;
		
        fsad += -fourpi*(*(Rho_Disk_Const_d+i))*f_prime*r*tlncoshz;
        fzad += -fourpi*(*(Rho_Disk_Const_d+i))*
                (f_prime*z*tlncoshz + f/(*(Z_Disk+i))*tanh(z/(*(Z_Disk+i))));
    }
}
        
__host__ __device__     
double AppDiskPotCUDA(double &r, double &z)
{
    //return 0;
    double appdiskpot = 0;
    double radius = sqrt(r*r+z*z);
    
    //Now add up the appdiskpot contributions from each sech^2 component
    //profile[3*1] is the density of the i-th component
    for (int i = 0; i < disk_d; ++i)
    {
        double Sigma_Profile[2], Rho_Profile[2];

        DiskProfileCUDA(radius, z, i, Sigma_Profile, Rho_Profile);
        
        double trunc_fac = GetTruncCUDA(radius, i);
        double zz = fabs(z/(*(Z_Disk+i)));
        appdiskpot += -fourpi*(*(Rho_Disk_Const_d+i))*Sigma_Profile[0]*(*(Z_Disk+i))*
                      (*(Z_Disk+i))*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
    }
    
    return appdiskpot;
}

__host__ __device__     
double AppDiskDensCUDA(double &r, double &z)
{
    //return 0;
    double appdiskdens = 0;
    double radius = sqrt(r*r+z*z);
    
    for (int i = 0; i < disk_d; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrimeCUDA(radius, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[3], Rho_Profile[3];

        DiskProfile2PrimeCUDA(radius, z, i, Sigma_Profile, Rho_Profile);

        double f = Sigma_Profile[0]*(*(Z_Disk+i))*(*(Z_Disk+i))*trunc_fac;
        double f1r = (*(Z_Disk+i))*(*(Z_Disk+i))*
                     (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/radius;
        double f2 = (*(Z_Disk+i))*(*(Z_Disk+i))*(Sigma_Profile[2]*trunc_fac +
                    2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		            (radius-(*(Out_Disk+i)))/(*(Dr_Trunc+i))/(*(Dr_Trunc+i)));

        if (radius==0)
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
    
    return appdiskdens;
}

__host__ __device__     
double DiskDensfICUDA(double r, double z, int i, thrust::device_ptr<double> Radius_CUDA, 
                                                 thrust::device_ptr<double> A_Pot_CUDA)
{
    double psi = PotCUDA(r, z, Radius_CUDA, A_Pot_CUDA);
    
    if(fabs(z/(*(Z_Disk+i)))>30)
        return 0;

    double trunc_fac = GetTruncCUDA(r, i);

    if (trunc_fac==0)
        return 0;

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
        double psizh = PotCUDA(r, 3.0*(*(Z_Disk+i)),  
                               Radius_CUDA, A_Pot_CUDA);
        double psi00 = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
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

//     double zfac = (*(Z_Disk+i))*(1-0.95*exp(-r*r/25));//+0.00001;
//     con = 1.0/cosh(fabs(z/zfac))/cosh(fabs(z/zfac));

	double Sigma_Profile[2], Rho_Profile[2];

    DiskProfileCUDA(r, z, i, Sigma_Profile, Rho_Profile);
    double disk_dens = (*(Rho_Disk_Const_d+i))*Sigma_Profile[0]*con*trunc_fac;
    
    return disk_dens;
}

