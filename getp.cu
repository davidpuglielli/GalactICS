struct GasDiskPolytropeFunctor
{
    int j;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> Polytrope_Const_CUDA;
    thrust::device_ptr<double> GasDensity_Const_CUDA;
    
    GasDiskPolytropeFunctor(int _j,
                            thrust::device_ptr<double> _Radius_CUDA,
                            thrust::device_ptr<double> _A_Pot_CUDA,
                            thrust::device_ptr<double> _Polytrope_Const_CUDA,
                            thrust::device_ptr<double> _GasDensity_Const_CUDA):
                            j(_j),
                            Radius_CUDA(_Radius_CUDA),
                            A_Pot_CUDA(_A_Pot_CUDA),
                            Polytrope_Const_CUDA(_Polytrope_Const_CUDA),
                            GasDensity_Const_CUDA(_GasDensity_Const_CUDA){}
    
    __host__ __device__ //__noinline__
    float operator()(double &radius) const
    {
        double poly_const = *(Sigma_0_Gas+j)*(*(Sigma_0_Gas+j));
        
//         if (radius==20||radius==10||radius==0)
//             printf("Did %d calculation      %f %f      %f\n", radius*100, radius,
//                    poly_const, PotCUDA(radius, 0, Radius_CUDA, A_Pot_CUDA));  
        
        return poly_const;
    }
};
    
struct GasDiskDensityFunctor
{
    int j;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> Polytrope_Const_CUDA;
    thrust::device_ptr<double> GasDensity_Const_CUDA;
    
    GasDiskDensityFunctor(int _j,
                          thrust::device_ptr<double> _Radius_CUDA,
                          thrust::device_ptr<double> _A_Pot_CUDA,
                          thrust::device_ptr<double> _Polytrope_Const_CUDA,
                          thrust::device_ptr<double> _GasDensity_Const_CUDA):
                          j(_j),
                          Radius_CUDA(_Radius_CUDA),
                          A_Pot_CUDA(_A_Pot_CUDA),
                          Polytrope_Const_CUDA(_Polytrope_Const_CUDA),
                          GasDensity_Const_CUDA(_GasDensity_Const_CUDA){}
    
    __host__ __device__ //__noinline__
    float operator()(double &radius) const
    {
        double gasdens_const = IterateGasDensCUDA(radius, j, Radius_CUDA, A_Pot_CUDA,
                                                  Polytrope_Const_CUDA, GasDensity_Const_CUDA);
        
        //if (radius%10<1e-6)
        //if (radius==20||radius==10||radius==0)
        //    printf("Did %d calculation      %f %f      %f\n", radius*100, radius,
        //           gasdens_const, PotCUDA(radius, 0, Radius_CUDA, A_Pot_CUDA));  
        
        return gasdens_const;
    }
};
    
__host__ __device__     
double GasDiskSurfaceDensfICUDA(double r, int j, 
                                thrust::device_ptr<double> Radius_CUDA, 
                                thrust::device_ptr<double> A_Pot_CUDA, 
                                thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    double dz = 0.025, zmax = 5;
    int integ_pts = int(zmax/dz);  
    double psi = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
    
    double plane_density = GasDiskDensfINoTruncCUDA(r, 0, j, Radius_CUDA, A_Pot_CUDA,
                                                    Polytrope_Const_CUDA, GasDensity_Const_CUDA)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    
    for (int i = 1; i < integ_pts; ++i)
    {
        //psi = PotCUDA(r, i*dz, Radius_CUDA, A_Pot_CUDA);
        //if(poly_const<10)
        //  cout << "iterate1 " << r << " " << i*dz << " " 
        //       << density << " " << psi << endl;
        density += GasDiskDensfINoTruncCUDA(r, i*dz, j, Radius_CUDA, A_Pot_CUDA,
                                            Polytrope_Const_CUDA, GasDensity_Const_CUDA)*dz;
    }
    
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << poly_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
__host__ __device__     
double GasDiskSurfaceDensfICUDA2(double r, int j, double dens_const,
                                 thrust::device_ptr<double> Radius_CUDA, 
                                 thrust::device_ptr<double> A_Pot_CUDA, 
                                 thrust::device_ptr<double> Polytrope_Const_CUDA, 
                                 thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    double dz = 0.025, zmax = 5;
    int integ_pts = int(zmax/dz);  
    double psi = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA);
    
    double plane_density = GasDiskDensfINoTruncCUDA2(r, 0, j, dens_const, Radius_CUDA, A_Pot_CUDA,
                                                     Polytrope_Const_CUDA, GasDensity_Const_CUDA)*dz;
    double density = 0;
    //cout << "iterate2 " << r << " " << 0 << " " << plane_density << " " << psi << endl;
    //printf("surfden %f %f %f\n", r, psi, plane_density);
    
    for (int i = 1; i < integ_pts; ++i)
    {
        //psi = PotCUDA(r, i*dz, Radius_CUDA, A_Pot_CUDA);
        //if(poly_const<10)
        //  cout << "iterate1 " << r << " " << i*dz << " " 
        //       << density << " " << psi << endl;
        density += GasDiskDensfINoTruncCUDA2(r, i*dz, j, dens_const, Radius_CUDA, A_Pot_CUDA,
                                            Polytrope_Const_CUDA, GasDensity_Const_CUDA)*dz;
    }
    
    //if(r<2&&r>1)
    //    printf("surfden %f %f %f\n", r, psi, density);
    density = 2*density+plane_density;
    
    //cout << "iterate " << r << " " << poly_const << " " << density << endl;
    //" " << target_density << endl;
    
    return density;
}
        
__host__ __device__     
double IterateGasDensCUDA(double r, int j,
                          thrust::device_ptr<double> Radius_CUDA, 
                          thrust::device_ptr<double> A_Pot_CUDA, 
                          thrust::device_ptr<double> Polytrope_Const_CUDA, 
                          thrust::device_ptr<double> GasDensity_Const_CUDA)
{
    //cout << "iterate " << endl;
    //if(r<100)
        //printf("iterate %f\n", r);
    
    double Sigma_Profile[2], Rho_Profile[2], zz = 0;
		
    GasDiskProfileCUDA(r, zz, j, Sigma_Profile, Rho_Profile);

    float target_density = *(GasDisk_Const_d+j)*Sigma_Profile[0]; 
       
    //cout << "iterate " << endl;
    //Find brackets for the polytropic constant. Here K is a function of density
    //so the following algorithm is a transcription of the root bisection
    //found in gendf.cpp. First, the find brackets part
    float dens_lower, dens_upper;
    dens_lower = dens_upper = 0.1;
    
    //cout << "iterate " << endl;
    //printf("iterate2\n");
    float rho_upper = GasDiskSurfaceDensfICUDA2(r, j, dens_upper, Radius_CUDA, 
                                                 A_Pot_CUDA, Polytrope_Const_CUDA, 
                                                 GasDensity_Const_CUDA);
    //cout << "iterate " << endl;
    float rho_lower = GasDiskSurfaceDensfICUDA2(r, j, dens_lower, Radius_CUDA, 
                                                 A_Pot_CUDA, Polytrope_Const_CUDA, 
                                                 GasDensity_Const_CUDA);
    
    //cout << "iterate " << endl;
    //printf("iterate3\n");
    while (rho_lower >= target_density)
    {
        dens_lower *= 0.1;
        rho_lower = GasDiskSurfaceDensfICUDA2(r, j, dens_lower, Radius_CUDA, 
                                              A_Pot_CUDA, Polytrope_Const_CUDA, 
                                              GasDensity_Const_CUDA);
        
        //printf("iterate3 %f\n", dens_lower);
        //printf("iterate4 %f %f\n", r, dens_upper);
        //printf("iterate %f\n", r);
        if (rho_lower == 0)
            return 0;
    }
    
    while (rho_upper <= target_density)
    {
        dens_upper *= 10;
        rho_upper = GasDiskSurfaceDensfICUDA2(r, j, dens_upper, Radius_CUDA, 
                                              A_Pot_CUDA, Polytrope_Const_CUDA, 
                                              GasDensity_Const_CUDA);
        
        //printf("iterate4 %f %f\n", r, dens_upper);
        //printf("iterate1 %f\n", r);
        if (rho_upper == 0)
            return 0;
    }
    
    //printf("iterate4\n");
    //Now the root bisection part
    rho_upper = GasDiskSurfaceDensfICUDA2(r, j, dens_upper, Radius_CUDA, 
                                          A_Pot_CUDA, Polytrope_Const_CUDA, 
                                          GasDensity_Const_CUDA) - target_density;
    rho_lower = GasDiskSurfaceDensfICUDA2(r, j, dens_lower, Radius_CUDA, 
                                          A_Pot_CUDA, Polytrope_Const_CUDA, 
                                          GasDensity_Const_CUDA) - target_density;
    
    //printf("iterate1\n");
    float dens_midpoint, rho_midpoint, poly_gap;
    int iters = 0;
    
    if (rho_upper*rho_lower > 0)
    {
        //cerr << "Bisection endpoints do not bracket root. Exiting..." << endl;
        //cerr << rho_upper << " " << rho_lower << " " << dens_lower << " " 
        //     << dens_upper << endl;
        exit(1);
    }
    
    do
    {
        ++iters;//cout<<"here "<<endl;
        
        dens_midpoint = (dens_upper+dens_lower)*0.5;
        
        rho_midpoint = GasDiskSurfaceDensfICUDA2(r, j, dens_midpoint, Radius_CUDA, 
                                                 A_Pot_CUDA, Polytrope_Const_CUDA, 
                                                 GasDensity_Const_CUDA) - target_density;
        
        if (rho_lower*rho_midpoint < 0)
        {
            dens_upper = dens_midpoint;
            rho_upper = rho_midpoint;
        }
        else if (rho_upper*rho_midpoint < 0)
        {
            dens_lower = dens_midpoint;
            rho_lower = rho_midpoint;
        }
        else if (rho_midpoint==0)
        {
            //cout << "\nFinished root bisection\n" << endl;
            return dens_midpoint;
        }
        else
        {
            //cerr << "Root bisection for dens failed! exiting..." << endl;
            //cerr << rho_lower << " " << rho_upper << " " << rho_midpoint << " "
            //     << dens_lower << " " << dens_upper << " " 
            //     << dens_midpoint << " " << target_density << endl;
            exit(1);
        }
        
        poly_gap = dens_upper-dens_lower;
    }
    while (iters < 200 && poly_gap > 0.01*dens_midpoint);
    
    //cout << "                                    " << iters << " " 
    //     << target_density << " " << dens_midpoint << endl;
    //if(r<10)printf("       %f %d %.12f    %e %e\n", r, iters, target_density, dens_midpoint, poly_gap);
    
    return dens_midpoint;
}
