__host__ __device__     
void GasDiskProfileCUDA(const double radius, const double z, const int i, 
                        double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0;//, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_GasDisk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        //ddensity = -density / G.(*(R_GasDisk+i));
        //dddensity = -ddensity / G.(*(R_GasDisk+i));
        //rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
        //              cosh(fabs(z/(*(Z_GasDisk+i)))));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        //ddensity = -density / G.(*(R_GasDisk+i));
        //dddensity = -ddensity / G.(*(R_GasDisk+i));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
                      cosh(fabs(z/(*(Z_GasDisk+i)))));
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
void GasDiskProfilePrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0;//, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_GasDisk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        ddensity = -density / (*(R_GasDisk+i));
        //dddensity = -ddensity / (*(R_GasDisk+i));
        //rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
        //              cosh(fabs(z/(*(Z_GasDisk+i)))));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        ddensity = -density / (*(R_GasDisk+i));
        //dddensity = -ddensity / (*(R_GasDisk+i));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
                      cosh(fabs(z/(*(Z_GasDisk+i)))));
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
void GasDiskProfile2PrimeCUDA(const double &radius, const double &z, const int &i, 
                           double *Sigma_Profile, double *Rho_Profile)
{
    double density = 0, ddensity = 0, dddensity = 0;
    double rho_density = 0;
    double zfac = (*(Z_GasDisk+i));//*(1-0.95*exp(-radius*radius/25));//+0.075;
    
    //Pure exponential
    if (i==0)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        ddensity = -density / (*(R_GasDisk+i));
        dddensity = -ddensity / (*(R_GasDisk+i));
        //rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
        //              cosh(fabs(z/(*(Z_GasDisk+i)))));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/zfac))*cosh(fabs(z/zfac)));
    }
    else if (i==1)
    {
        density = exp(-radius/(*(R_GasDisk+i)));
        ddensity = -density / (*(R_GasDisk+i));
        dddensity = -ddensity / (*(R_GasDisk+i));
        rho_density = exp(-radius/(*(R_GasDisk+i)))/(cosh(fabs(z/(*(Z_GasDisk+i))))*
                      cosh(fabs(z/(*(Z_GasDisk+i)))));
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
void AppGasDiskForceCUDA(double &r, double &z, double &fsad, double &fzad)
{
    double s = sqrt(r*r + z*z);
    fsad = fzad = 0;
    return;
        
    for (int i = 0; i < disk_d; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrimeCUDA(s, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[2], Rho_Profile[2];

        GasDiskProfilePrimeCUDA(s, z, i, Sigma_Profile, Rho_Profile);
    
        double zz = fabs(z/(*(Z_GasDisk+i)));
        double e2zz = exp(-2*zz);
        double tlncoshz = zz+log(0.5*(1.0 + e2zz));
        
		double f = Sigma_Profile[0]*trunc_fac;
		double f_prime = Sigma_Profile[0]*trunc_facprime + 
                         Sigma_Profile[1]*trunc_fac;
		
        fsad += -fourpi*(*(Rho_GasDisk_Const_d+i))*f_prime*r*tlncoshz;
        fzad += -fourpi*(*(Rho_GasDisk_Const_d+i))*
                (f_prime*z*tlncoshz + f/(*(Z_GasDisk+i))*tanh(z/(*(Z_GasDisk+i))));
    }
}
        
__host__ __device__     
double AppGasDiskPotCUDA(double &r, double &z)
{
    return 0;
    double appdiskpot = 0;
    double radius = sqrt(r*r+z*z);
    
    //Now add up the appdiskpot contributions from each sech^2 component
    //profile[3*1] is the density of the i-th component
    for (int i = 0; i < disk_d; ++i)
    {
        double Sigma_Profile[2], Rho_Profile[2];

        GasDiskProfileCUDA(radius, z, i, Sigma_Profile, Rho_Profile);
        
        double trunc_fac = GetTruncCUDA(radius, i);
        double zz = fabs(z/(*(Z_GasDisk+i)));
        appdiskpot += -fourpi*(*(Rho_GasDisk_Const_d+i))*Sigma_Profile[0]*(*(Z_GasDisk+i))*
                      (*(Z_GasDisk+i))*(zz + log(0.5 + 0.5*exp(-2*zz))) * trunc_fac;
    }
    
    return appdiskpot;
}

__host__ __device__     
double AppGasDiskDensCUDA(double &r, double &z)
{
    return 0;
    double appdiskdens = 0;
    double radius = sqrt(r*r+z*z);
    
    for (int i = 0; i < disk_d; ++i)
    {
        double trunc_fac, trunc_facprime;

        GetTruncPrimeCUDA(radius, i, trunc_fac, trunc_facprime);

        double Sigma_Profile[3], Rho_Profile[3];

        GasDiskProfile2PrimeCUDA(radius, z, i, Sigma_Profile, Rho_Profile);

        double f = Sigma_Profile[0]*(*(Z_GasDisk+i))*(*(Z_GasDisk+i))*trunc_fac;
        double f1r = (*(Z_GasDisk+i))*(*(Z_GasDisk+i))*
                     (Sigma_Profile[1]*trunc_fac + Sigma_Profile[0]*trunc_facprime)/radius;
        double f2 = (*(Z_GasDisk+i))*(*(Z_GasDisk+i))*(Sigma_Profile[2]*trunc_fac +
                    2*Sigma_Profile[1]*trunc_facprime - Sigma_Profile[0]*trunc_facprime*
		            (radius-(*(Out_GasDisk+i)))/(*(Dr_Trunc+i))/(*(Dr_Trunc+i)));

        if (radius==0)
        {
            f1r = f2 = 0;
        }        

        double zz = fabs(z/(*(Z_GasDisk+i)));
        double ezz = exp(-zz), e2zz = exp(-2*zz);
        double tlncosh = zz+log(0.5*(1.0 + e2zz));
        double tztanh = zz*(1-e2zz)/(1+e2zz);
        double tsech2 = (2*ezz/(1+e2zz))*(2*ezz/(1+e2zz));
        double total = f2*tlncosh + 2*f1r*(tztanh+tlncosh) + 
                       f*tsech2/(*(Z_GasDisk+i))/(*(Z_GasDisk+i));
        appdiskdens += (*(Rho_GasDisk_Const_d+i))*total;
    }
    
    return appdiskdens;
}

__host__ __device__     
double GasDiskDensfICUDA(double r, double z, int i, 
                         thrust::device_ptr<double> Radius_CUDA, 
                         thrust::device_ptr<double> A_Pot_CUDA, 
                         thrust::device_ptr<double> Polytrope_Const_CUDA, 
                         thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    double psi = PotCUDA(r, z, Radius_CUDA, A_Pot_CUDA);
    double psi_0 = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
    psi -= psi_0;
    
    double trunc_fac = GetTruncGasCUDA(r, i);

       // printf("gasdens    %f %f %f %f\n", r, z, psi, trunc_fac);
    if (trunc_fac==0)
    {
        return 0;
    }
    else if (*(Gamma+i) == 1)
    {
        if(fabs(z/Z_GasDisk[i])>30)
            return 0;

        double trunc_fac = GetTruncGasCUDA(r, i);

        if (trunc_fac==0)
            return 0;

        double poly_const = GetPolyConstCUDA(r, i, Radius_CUDA, Polytrope_Const_CUDA);
        double dens_const = GetDensConstCUDA(r, i, Radius_CUDA, GasDensity_Const_CUDA);

	    //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density, rho_0 = dens_const;//*Sigma_Profile[0];

        density = rho_0*exp(psi/poly_const)*trunc_fac;//*Sigma_Profile[0];

        //if(r<10)printf("gasdens    %f %f %f %f\n", r, z, psi, density);
        
        if (density < 0)
        {
            return 0;
        }

        return density;
    }
    else
    {
        return 0;
    }
}

__host__ __device__     
double GasDiskDensfINoTruncCUDA(double r, double z, int i, 
                                thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> A_Pot_CUDA, 
                                thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    double psi = PotCUDA(r, z, Radius_CUDA, A_Pot_CUDA);
    double psi_0 = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
    psi -= psi_0;
    
    double trunc_fac = GetTruncGasCUDA(r, i);

    if (trunc_fac==0)
    {
        return 0;
    }
    if (*(Gamma+i) == 1)
    {
        if(fabs(z/Z_GasDisk[i])>30)
            return 0;

        double trunc_fac = GetTruncGasCUDA(r, i);

        if (trunc_fac==0)
            return 0;

        double poly_const = GetPolyConstCUDA(r, i, Radius_CUDA, Polytrope_Const_CUDA);
        double dens_const = GetDensConstCUDA(r, i, Radius_CUDA, GasDensity_Const_CUDA);

	    //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density, rho_0 = dens_const*trunc_fac;//*Sigma_Profile[0];

        density = rho_0*exp(psi/poly_const);

        if (density < 0)
        {
            return 0;
        }

        return density;
    }
    else
    {
        return 0;
    }
}

__host__ __device__     
double GasDiskDensfINoTruncCUDA2(double r, double z, int i, double dens_const,
                                 thrust::device_ptr<double> Radius_CUDA, 
                                 thrust::device_ptr<double> A_Pot_CUDA, 
                                 thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                 thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    double psi = PotCUDA(r, z, Radius_CUDA, A_Pot_CUDA);
    double psi_0 = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
    psi -= psi_0;
    
    //if (r<1) printf("den  %f %f %f %f\n", r, z, psi, *(Gamma+i));
    
    double trunc_fac = GetTruncGasCUDA(r, i);

    if (trunc_fac==0)
    {
        return 0;
    }
    if (*(Gamma+i) == 1)
    {
        if(fabs(z/Z_GasDisk[i])>30)
            return 0;

        double trunc_fac = GetTruncGasCUDA(r, i);

        if (trunc_fac==0)
            return 0;

        double poly_const = GetPolyConstCUDA(r, i, Radius_CUDA, Polytrope_Const_CUDA);
        //double dens_const = GetDensConstCUDA(r, i, Radius_CUDA, GasDensity_Const_CUDA);

	    //double Sigma_Profile[2], Rho_Profile[2];

        //GasDiskProfile(r, z, i, Sigma_Profile, Rho_Profile);

        double density, rho_0 = dens_const;//*Sigma_Profile[0];

        density = rho_0*exp(psi/poly_const)*trunc_fac;
        //if (r>1&&r<2) printf("den  %f %f %f   %f %f\n", r, z, psi, dens_const, poly_const);

        if (density < 0)
        {
            return 0;
        }

        return density;
    }
    else
    {
        return 0;
    }
}

__host__ __device__     
double GetDensConstCUDA(double r, int j, thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> GasDensity_Const_CUDA)
{
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
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_d-1;
    }
    
    double r1 = Radius_CUDA[ihi-1];
    double r2 = Radius_CUDA[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*GasDensity_Const_CUDA[j*nr_d+ihi] + tm1*GasDensity_Const_CUDA[j*nr_d+ihi-1];
}
    
__host__ __device__     
double GetPolyConstCUDA(double r, int j, thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> Polytrope_Const_CUDA)
{
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
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_d-1;
    }
    
    double r1 = Radius_CUDA[ihi-1];
    double r2 = Radius_CUDA[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*Polytrope_Const_CUDA[j*nr_d+ihi] + tm1*Polytrope_Const_CUDA[j*nr_d+ihi-1];
}
    
// __host__ __device__     
// double GetGasHeightCUDA(double r, int j, thrust::device_ptr<double> Radius_CUDA,
//                         thrust::device_ptr<double> GasDensity_Height_CUDA)
// {
//     int ihi = ceil(r/dr_d);
//     
//     if(r < dr_d)
//     {
//         ihi = 1;
//     }
//     else if (ihi < 1)
//     {
//     }
//     else if (ihi > nr_d-1)
//     {
//         //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
//         //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
//         ihi = nr_d-1;
//     }
//     
//     double r1 = Radius_CUDA[ihi-1];
//     double r2 = Radius_CUDA[ihi];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     return t*GasDensity_Height_CUDA[j*nr_d+ihi] + tm1*GasDensity_Height_CUDA[j*nr_d+ihi-1];
// }
//     
//     
//     
