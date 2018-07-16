struct SurfaceBrightnessFunctor
{
    const double error_factor;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> F_R_CUDA;
    thrust::device_ptr<double> DF_Sersic_CUDA;
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;

    SurfaceBrightnessFunctor(double _error_factor,
                             thrust::device_ptr<double> _Radius_CUDA,
                             thrust::device_ptr<double> _A_Pot_CUDA,
                             thrust::device_ptr<double> _F_R_CUDA,
                             thrust::device_ptr<double> _DF_Sersic_CUDA,
                             thrust::device_ptr<double> _Dens_Psi_Bulge_CUDA): 
                             error_factor(_error_factor),
                             Radius_CUDA(_Radius_CUDA),
                             A_Pot_CUDA(_A_Pot_CUDA),
                             F_R_CUDA(_F_R_CUDA),
                             DF_Sersic_CUDA(_DF_Sersic_CUDA),
                             Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge_CUDA){}
    //StarVelocityFunctor(){}

    //SurfaceBrightnessFunctor(double _error_factor): error_factor(_error_factor){}
    //SurfaceBrightnessFunctor(){}
    
    template <typename Tuple>
    __host__ __device__ //__noinline__ 
    float operator()(const Tuple data_point) const
    {
        double radius       = thrust::get<0>(data_point);
        double data         = thrust::get<1>(data_point);
        double error        = thrust::get<2>(data_point);
        double ellipticity  = thrust::get<3>(data_point);
        
        double integral_bulge = 0;
        double *integral_disk = new double[disk_d], *disk_light = new double[disk_d];
        double kpc_to_arcsec = atan(1.0/distance_d)*3600*180/PI;//arcsec per kpc
        
        error = sqrt(error*error+error_factor*error_factor);
                
        double point_inclination = acos(1.0-ellipticity);
        
        double d_phi = PI/2/99;
        
        //Get force
        //double f_rad, f_z;
        //ForceCUDA(radius, 0, f_rad, f_z, Radius_CUDA, A_Pot_CUDA, F_R_CUDA);
        
        for (int j = 0; j < disk_d; ++j)
        {
            integral_disk[j] = 0;
            disk_light[j] = 0;
        }
    
        for (int j = 0; j < 100; ++j)
        {
            double yp = radius*cos(j*d_phi);
            double xp = radius*sin(j*d_phi)*cos(point_inclination);
            //cout << "step " << i << " " << j << " " << xp << " " << yp << endl;
            
            for (int k = 0; k < disk_d; ++k)
            {
                double disk_surf_den = DiskSurfaceDensityCUDA(xp, yp, k, Radius_CUDA, A_Pot_CUDA);
                integral_disk[k] += disk_surf_den;
            }
            
            double bulge_surf_den = BulgeSurfaceDensityCUDA(xp, yp, Radius_CUDA,
                                                            A_Pot_CUDA, Dens_Psi_Bulge_CUDA);
            
            if (j == 0 || j == 99)
            {
                bulge_surf_den *= 0.5;
            }
            
            integral_bulge += bulge_surf_den;
        }
        
        //printf("in functor3\n");
        double bulge_light = 0, disk_light_tot = 0;
        
        for (int j = 0; j < disk_d; ++j)
        {
            integral_disk[j] /= 100;
            disk_light[j] += integral_disk[j]/(*(ML_Disk+j))*m_scale/kpc_to_arcsec/kpc_to_arcsec;
            disk_light_tot += disk_light[j];
        }
        
        integral_bulge /= 100;
        bulge_light = integral_bulge/ml_bulge_d*m_scale/kpc_to_arcsec/kpc_to_arcsec;
        
        //the 1e-6 is for conv to pc in log arg
        double disk_magnitude_tot = mag_sun-5-2.5*log10(disk_light_tot/distance_d/distance_d*1e-6);
        double bulge_magnitude = mag_sun-5-2.5*log10(bulge_light/distance_d/distance_d*1e-6);
        double total_magnitude = mag_sun-5-2.5*log10((disk_light_tot+bulge_light)/distance_d/distance_d*1e-6);
        double *disk_magnitude = new double[disk_d];
        
        for (int j = 0; j < disk_d; ++j)
        {
            disk_magnitude[j] = mag_sun-5-2.5*log10(disk_light[j]/distance_d/distance_d*1e-6);
        }
        
        //printf("in functor4 %f %f %f\n", total_magnitude, data, error);
        double chi_square = (total_magnitude-data)*(total_magnitude-data)/error/error + 2*log(error); 
        
        //printf("in functor5\n");
        return chi_square;       
    }
};

struct StarVelocityFunctor
{
    const double error_factor;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> Rad_Spline_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> F_R_CUDA;
    thrust::device_ptr<double> Omega_CUDA;
    thrust::device_ptr<double> Omega2_CUDA;
    thrust::device_ptr<double> A_K_CUDA;
    thrust::device_ptr<double> A_K2_CUDA;
    thrust::device_ptr<double> DF_Sersic_CUDA;
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
    thrust::device_ptr<double> Am_Tab_CUDA;
    thrust::device_ptr<double> R_Tab_CUDA;
    thrust::device_ptr<double> R_Tab2_CUDA;
    thrust::device_ptr<double> R_Tab2_Zero_CUDA;
    thrust::device_ptr<double> FD_Rat_CUDA;
    thrust::device_ptr<double> D_Rat2_CUDA;
    thrust::device_ptr<double> FSZ_Rat_CUDA;
    thrust::device_ptr<double> SZ_Rat2_CUDA;

    StarVelocityFunctor(double _error_factor,
                        thrust::device_ptr<double> _Radius_CUDA,
                        thrust::device_ptr<double> _Rad_Spline_CUDA,
                        thrust::device_ptr<double> _A_Pot_CUDA,
                        thrust::device_ptr<double> _F_R_CUDA,
                        thrust::device_ptr<double> _Omega_CUDA,
                        thrust::device_ptr<double> _Omega2_CUDA,
                        thrust::device_ptr<double> _A_K_CUDA,
                        thrust::device_ptr<double> _A_K2_CUDA,
                        thrust::device_ptr<double> _DF_Sersic_CUDA,
                        thrust::device_ptr<double> _Dens_Psi_Bulge_CUDA,
                        thrust::device_ptr<double> _Am_Tab_CUDA,
                        thrust::device_ptr<double> _R_Tab_CUDA,
                        thrust::device_ptr<double> _R_Tab2_CUDA,
                        thrust::device_ptr<double> _R_Tab2_Zero_CUDA,
                        thrust::device_ptr<double> _FD_Rat_CUDA,
                        thrust::device_ptr<double> _D_Rat2_CUDA,
                        thrust::device_ptr<double> _FSZ_Rat_CUDA,
                        thrust::device_ptr<double> _SZ_Rat2_CUDA):
                        error_factor(_error_factor),
                        Radius_CUDA(_Radius_CUDA),
                        Rad_Spline_CUDA(_Rad_Spline_CUDA),
                        A_Pot_CUDA(_A_Pot_CUDA),
                        F_R_CUDA(_F_R_CUDA),
                        Omega_CUDA(_Omega_CUDA),
                        Omega2_CUDA(_Omega2_CUDA),
                        A_K_CUDA(_A_K_CUDA),
                        A_K2_CUDA(_A_K2_CUDA),
                        DF_Sersic_CUDA(_DF_Sersic_CUDA),
                        Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge_CUDA),
                        Am_Tab_CUDA(_Am_Tab_CUDA),
                        R_Tab_CUDA(_R_Tab_CUDA),
                        R_Tab2_CUDA(_R_Tab2_CUDA),
                        R_Tab2_Zero_CUDA(_R_Tab2_Zero_CUDA),
                        FD_Rat_CUDA(_FD_Rat_CUDA),
                        D_Rat2_CUDA(_D_Rat2_CUDA),
                        FSZ_Rat_CUDA(_FSZ_Rat_CUDA),
                        SZ_Rat2_CUDA(_SZ_Rat2_CUDA){}
    //StarVelocityFunctor(){}
    
    //__host__ __device__ 
    //float operator()(tuple<double, double, double> data_point) const
    template <typename Tuple>
    __host__ __device__ //__noinline__ 
    float operator()(const Tuple data_point) const
    {
        double radius = thrust::get<0>(data_point);
        double data   = thrust::get<1>(data_point);
        double error  = thrust::get<2>(data_point);
        
        error = sqrt(error*error+error_factor*error_factor);
        
        //printf("functor1\n");
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        double x = 0;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < disk_d; ++j)
        {
            disk_surface_density1 = DiskSurfaceDensityCUDA(x, radius, j, Radius_CUDA, A_Pot_CUDA);
            disk_surface_density += disk_surface_density1;
            surface_light += disk_surface_density/(*(ML_Disk+j));
            
            DiskDispersionCUDA(x, radius, wbar_disk1, disk_dispersion1, j, Radius_CUDA, Rad_Spline_CUDA,
                               A_Pot_CUDA, F_R_CUDA, Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
                               A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA,
                               R_Tab2_Zero_CUDA, FD_Rat_CUDA,
                               D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
            wbar_disk += wbar_disk1/(*(ML_Disk+j))*disk_surface_density1;
        } 
        
        //Add on the bulge contribution
        double bulge_surface_density = 0;
        double wbar_bulge = 0, bulge_dispersion;
        
        if (bulge_flag_d)
        {
            bulge_surface_density = BulgeSurfaceDensityCUDA(x, radius, 
                                    Radius_CUDA, A_Pot_CUDA, Dens_Psi_Bulge_CUDA);
        
            surface_light += bulge_surface_density/ml_bulge_d;
                
            BulgeDispersionCUDA(x, radius, wbar_bulge, bulge_dispersion, 
                                Radius_CUDA, A_Pot_CUDA, Dens_Psi_Bulge_CUDA, DF_Sersic_CUDA);
        }
        
        wbar = (wbar_disk + wbar_bulge/ml_bulge_d*bulge_surface_density)/surface_light;
        
        double chi_square = (wbar-data)*(wbar-data)/error/error + 2*log(error); 
        
        return chi_square;       
    }
};

struct StarDispersionFunctor
{
    const double error_factor;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> Rad_Spline_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> F_R_CUDA;
    thrust::device_ptr<double> Omega_CUDA;
    thrust::device_ptr<double> Omega2_CUDA;
    thrust::device_ptr<double> A_K_CUDA;
    thrust::device_ptr<double> A_K2_CUDA;
    thrust::device_ptr<double> DF_Sersic_CUDA;
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
    thrust::device_ptr<double> Am_Tab_CUDA;
    thrust::device_ptr<double> R_Tab_CUDA;
    thrust::device_ptr<double> R_Tab2_CUDA;
    thrust::device_ptr<double> R_Tab2_Zero_CUDA;
    thrust::device_ptr<double> FD_Rat_CUDA;
    thrust::device_ptr<double> D_Rat2_CUDA;
    thrust::device_ptr<double> FSZ_Rat_CUDA;
    thrust::device_ptr<double> SZ_Rat2_CUDA;

    StarDispersionFunctor(double _error_factor,
                          thrust::device_ptr<double> _Radius_CUDA,
                          thrust::device_ptr<double> _Rad_Spline_CUDA,
                          thrust::device_ptr<double> _A_Pot_CUDA,
                          thrust::device_ptr<double> _F_R_CUDA,
                          thrust::device_ptr<double> _Omega_CUDA,
                          thrust::device_ptr<double> _Omega2_CUDA,
                          thrust::device_ptr<double> _A_K_CUDA,
                          thrust::device_ptr<double> _A_K2_CUDA,
                          thrust::device_ptr<double> _DF_Sersic_CUDA,
                          thrust::device_ptr<double> _Dens_Psi_Bulge_CUDA,
                          thrust::device_ptr<double> _Am_Tab_CUDA,
                          thrust::device_ptr<double> _R_Tab_CUDA,
                          thrust::device_ptr<double> _R_Tab2_CUDA,
                          thrust::device_ptr<double> _R_Tab2_Zero_CUDA,
                          thrust::device_ptr<double> _FD_Rat_CUDA,
                          thrust::device_ptr<double> _D_Rat2_CUDA,
                          thrust::device_ptr<double> _FSZ_Rat_CUDA,
                          thrust::device_ptr<double> _SZ_Rat2_CUDA):
                          error_factor(_error_factor),
                          Radius_CUDA(_Radius_CUDA),
                          Rad_Spline_CUDA(_Rad_Spline_CUDA),
                          A_Pot_CUDA(_A_Pot_CUDA),
                          F_R_CUDA(_F_R_CUDA),
                          Omega_CUDA(_Omega_CUDA),
                          Omega2_CUDA(_Omega2_CUDA),
                          A_K_CUDA(_A_K_CUDA),
                          A_K2_CUDA(_A_K2_CUDA),
                          DF_Sersic_CUDA(_DF_Sersic_CUDA),
                          Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge_CUDA),
                          Am_Tab_CUDA(_Am_Tab_CUDA),
                          R_Tab_CUDA(_R_Tab_CUDA),
                          R_Tab2_CUDA(_R_Tab2_CUDA),
                          R_Tab2_Zero_CUDA(_R_Tab2_Zero_CUDA),
                          FD_Rat_CUDA(_FD_Rat_CUDA),
                          D_Rat2_CUDA(_D_Rat2_CUDA),
                          FSZ_Rat_CUDA(_FSZ_Rat_CUDA),
                          SZ_Rat2_CUDA(_SZ_Rat2_CUDA){}
    //StarDispersionFunctor(){}
    
    //__host__ __device__ 
    //float operator()(tuple<double, double, double> data_point) const
    template <typename Tuple>
    __host__ __device__ //__noinline__ 
    float operator()(const Tuple data_point) const
    {
        double radius = thrust::get<0>(data_point);
        double data   = thrust::get<1>(data_point);
        double error  = thrust::get<2>(data_point);
        
        error = sqrt(error*error+error_factor*error_factor);
        
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, disk_dispersion = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        double x = 0;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < disk_d; ++j)
        {
            disk_surface_density1 = DiskSurfaceDensityCUDA(x, radius, j, Radius_CUDA, A_Pot_CUDA);
            disk_surface_density += disk_surface_density1;
            surface_light += disk_surface_density/(*(ML_Disk+j));
            
            DiskDispersionCUDA(x, radius, wbar_disk1, disk_dispersion1, j, Radius_CUDA, Rad_Spline_CUDA,
                               A_Pot_CUDA, F_R_CUDA, Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
                               A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA,
                               R_Tab2_Zero_CUDA, FD_Rat_CUDA,
                               D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
            wbar_disk += wbar_disk1/(*(ML_Disk+j))*disk_surface_density1;
            disk_dispersion += disk_dispersion1/(*(ML_Disk+j))*disk_surface_density1;
        } 
        
        //Add on the bulge contribution
        double bulge_surface_density = 0;
        double wbar_bulge = 0, bulge_dispersion = 0, total_dispersion;
        double disk_surface_light = surface_light, bulge_surface_light = 0;
        
        if (bulge_flag_d)
        {
            bulge_surface_density = BulgeSurfaceDensityCUDA(x, radius, 
                                    Radius_CUDA, A_Pot_CUDA, Dens_Psi_Bulge_CUDA);
        
            bulge_surface_light = bulge_surface_density/ml_bulge_d;
            surface_light += bulge_surface_light;
                
            BulgeDispersionCUDA(x, radius, wbar_bulge, bulge_dispersion, 
                                Radius_CUDA, A_Pot_CUDA, Dens_Psi_Bulge_CUDA, DF_Sersic_CUDA);
        }
        
        wbar = (wbar_disk + wbar_bulge/ml_bulge_d*bulge_surface_density)/surface_light;
        
        total_dispersion = (disk_dispersion + bulge_dispersion/ml_bulge_d*bulge_surface_density)/
                           surface_light;
                           
        total_dispersion = sqrt(total_dispersion - wbar*wbar);
        
        disk_dispersion = sqrt(disk_dispersion/disk_surface_light - 
                               wbar_disk*wbar_disk/disk_surface_light/disk_surface_light);
        
        if (bulge_flag_d)
        {
            bulge_dispersion = sqrt(bulge_dispersion - wbar_bulge*wbar_bulge);
        }
        
        double chi_square = (total_dispersion-data)*(total_dispersion-data)/error/error + 2*log(error); 
        
        return chi_square;       
    }
};

struct CircVelocityFunctor
{
    const double error_factor;
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> A_Pot_CUDA;
    thrust::device_ptr<double> F_R_CUDA;

    CircVelocityFunctor(double _error_factor,
                        thrust::device_ptr<double> _Radius_CUDA,
                        thrust::device_ptr<double> _A_Pot_CUDA,
                        thrust::device_ptr<double> _F_R_CUDA): 
                        error_factor(_error_factor), 
                        Radius_CUDA(_Radius_CUDA),
                        A_Pot_CUDA(_A_Pot_CUDA), 
                        F_R_CUDA(_F_R_CUDA){}
    //StarVelocityFunctor(){}
    
    template <typename Tuple>
    __host__ __device__ //__noinline__ 
    float operator()(const Tuple data_point) const
    {
        double radius = thrust::get<0>(data_point);
        double data   = thrust::get<1>(data_point);
        double error  = thrust::get<2>(data_point);
        
        error = sqrt(error*error+error_factor*error_factor);
        
        //Get force
        double f_rad, f_z;
        ForceCUDA(radius, 0, f_rad, f_z, Radius_CUDA, A_Pot_CUDA, F_R_CUDA);
        
        double v_circ = sqrt(-radius*f_rad)*sin_inclination_d;
        
        double chi_square = (v_circ-data)*(v_circ-data)/error/error+2*log(error);
        
        return chi_square;
    }
};

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Chi Square functions /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

__host__ __device__     
void RotateCoordinatesCUDA(double &xp, double &yp, double &zp, 
                           double &x, double &y, double &z)
{
    x = xp*cos_inclination_d - zp*sin_inclination_d;
    y = yp;
    z = xp*sin_inclination_d + zp*cos_inclination_d;
}    
	
__host__ __device__     
void RotateCoordinatesBackCUDA(double &xp, double &yp, double &zp, 
                               double &x, double &y, double &z)
{
    x = xp*cos_inclination_d + zp*sin_inclination_d;
    y = yp;
    z = -xp*sin_inclination_d + zp*cos_inclination_d;
}    

__host__ __device__     
void RotateCoordinatesOffCentreCUDA(double &xp, double &yp, double &zp, 
                                    double &x, double &y, double &z)
{
    x = xp/cos_inclination_d - zp*sin_inclination_d;
    y = yp;
    z = zp*cos_inclination_d;
}    

__host__ __device__     
double DiskSurfaceDensityCUDA(double &xp, double &yp, int &j,
                              thrust::device_ptr<double> Radius_CUDA, 
                              thrust::device_ptr<double> A_Pot_CUDA)
{
    double x, y, z;
    double sum = 0, zdc = *(Z_Disk+j)/cos_inclination_d;
    double dt = 8*zdc/n_los_d;
    double zp;
    //double weights = 0;
    
    //RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
        //printf("SurfDens1\n");
    
    double r_cyl, disk_dens;
    //sum += 1.5*disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
    
    for (int i = 1; i < n_los_d; ++i)
    {
        zp = -4*zdc+i*dt;;
    
        RotateCoordinatesOffCentreCUDA(xp, yp, zp, x, y, z);
        //printf("SurfDens1 %f %f %f %f %f %f\n", xp, yp, zp, x, y, z);
    
        r_cyl = sqrt(x*x+y*y);
        disk_dens = DiskDensfICUDA(r_cyl, z, j, Radius_CUDA, A_Pot_CUDA);
        sum += disk_dens*dt;
        //weights += dt;
    }
    
    //Now do the i=1 and i=n_los-1 cases
    zp = -4*zdc;

        //printf("SurfDens2\n");
    RotateCoordinatesOffCentreCUDA(xp, yp, zp, x, y, z);

        //printf("SurfDens31\n");
    r_cyl = sqrt(x*x+y*y);
        //printf("SurfDens32\n");
    disk_dens = DiskDensfICUDA(r_cyl, z, j, Radius_CUDA, A_Pot_CUDA);
        //printf("SurfDens33\n");
    sum += 0.5*disk_dens*dt;
        //printf("SurfDens34\n");
    //weights += 0.5*cosh(zp/zdc)*cosh(zp/zdc);
    
    zp = 4*zdc;

    RotateCoordinatesOffCentreCUDA(xp, yp, zp, x, y, z);

        //printf("SurfDens35\n");
    r_cyl = sqrt(x*x+y*y);
        //printf("SurfDens36 %f %f %d\n", r_cyl, z, j);
    disk_dens = DiskDensfICUDA(r_cyl, z, j, Radius_CUDA, A_Pot_CUDA);
        //printf("SurfDens37\n");
    sum += 0.5*disk_dens*dt;
    //weights += 1.5*cosh(zp/zdc)*cosh(zp/zdc);
        //printf("SurfDens3\n");
    
    return sum  ;
}

__host__ __device__     
double BulgeSurfaceDensityCUDA(double &xp, double &yp, 
                               thrust::device_ptr<double> Radius_CUDA, 
                               thrust::device_ptr<double> A_Pot_CUDA,
                               thrust::device_ptr<double> Dens_Psi_Bulge_CUDA)
{
    double x, y, z;
    double sum = 0;
    double rp = xp*xp+yp*yp;
    double rlogmax = log(max(a_bulge_d*pow(10.0, n_sersic_d), r_max_d)), rlogmin;
    
    if (rp > 0)
    {
        rlogmin = log(rp/100);
    }
    else 
    {
        rlogmin = log(0.1*dr_d);
    }
    
    double drlog = (rlogmax-rlogmin)/n_los_d;
    double deltar = exp(rlogmin);
    
    double zp = 0;
    
    RotateCoordinatesCUDA(xp, yp, zp, x, y, z);
    
    double r_cyl = sqrt(x*x+y*y);
    double psi = PotCUDA(r_cyl, z, Radius_CUDA, A_Pot_CUDA);
    double bulge_dens = BulgeDensPsiCUDA(psi, Dens_Psi_Bulge_CUDA);

    sum += 0.5*deltar*bulge_dens;
    
    for (int i = 0; i < n_los_d; ++i)
    {
        double zlog = rlogmin + i*drlog;
        zp = exp(zlog);
        RotateCoordinatesOffCentreCUDA(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        psi = PotCUDA(r_cyl, z, Radius_CUDA, A_Pot_CUDA);
        bulge_dens = BulgeDensPsiCUDA(psi, Dens_Psi_Bulge_CUDA);
        sum += drlog*zp*bulge_dens;
    }
    
    return 2*sum;
}

__host__ __device__     
void PreDiskVelocitiesCUDA(double *pos, double &v_phi_max, double *v_cyl, 
                           double &f_max, double &v_circ, int &j, 
                           thrust::device_ptr<double> Radius_CUDA,
                           thrust::device_ptr<double> Rad_Spline_CUDA,
                           thrust::device_ptr<double> A_Pot_CUDA,
                           thrust::device_ptr<double> F_R_CUDA,
                           thrust::device_ptr<double> Omega_CUDA,
                           thrust::device_ptr<double> Omega2_CUDA,
                           thrust::device_ptr<double> A_K_CUDA,
                           thrust::device_ptr<double> A_K2_CUDA,
                           thrust::device_ptr<double> Am_Tab_CUDA,
                           thrust::device_ptr<double> R_Tab_CUDA,
                           thrust::device_ptr<double> R_Tab2_CUDA,
                           thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                           thrust::device_ptr<double> FD_Rat_CUDA,
                           thrust::device_ptr<double> D_Rat2_CUDA,
                           thrust::device_ptr<double> FSZ_Rat_CUDA,
                           thrust::device_ptr<double> SZ_Rat2_CUDA)
{
    double freq_omega, freq_kappa, x = pos[0], y = pos[1], z = pos[2];
    double r_cyl = sqrt(x*x+y*y);
    double fs, fz;
    GetOmegaKappaCUDA(r_cyl, freq_omega, freq_kappa, Radius_CUDA, Omega_CUDA,
                      Omega2_CUDA, A_K_CUDA, A_K2_CUDA);
    
    v_phi_max = r_cyl*freq_omega;
    v_cyl[0] = sqrt(SigR2CUDA(r_cyl, j));
    v_cyl[1] = 0.5*freq_kappa*v_cyl[0]/freq_omega;
    v_cyl[2] = sqrt(SigZ2CUDA(r_cyl, j, Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, 
                              FSZ_Rat_CUDA, SZ_Rat2_CUDA));
    
    double v_phi_max_old = v_phi_max;
    
    FindMaxCUDA(r_cyl, z, v_cyl[1], v_phi_max, f_max, j, Radius_CUDA, 
                Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
                A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA, 
                R_Tab2_Zero_CUDA, FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, 
                SZ_Rat2_CUDA);
    f_max *= 1.1;
    
    v_phi_max = v_phi_max_old;
    
    ForceCUDA(r_cyl, z, fs, fz, Radius_CUDA, A_Pot_CUDA, F_R_CUDA);
    v_circ = sqrt(-fs*r_cyl);
}
    
__host__ __device__     
void DiskVelocitiesCUDA(double *pos, double &v_phi_max, double *v_cyl, double &f_max,
                        double &v_circ, double &w_los, double &df, int &j, 
                        thrust::device_ptr<double> Radius_CUDA, 
                        thrust::device_ptr<double> Rad_Spline_CUDA,
                        thrust::device_ptr<double> A_Pot_CUDA,
                        thrust::device_ptr<double> Omega_CUDA,
                        thrust::device_ptr<double> Omega2_CUDA,
                        thrust::device_ptr<double> A_K_CUDA,
                        thrust::device_ptr<double> A_K2_CUDA,
                        thrust::device_ptr<double> Am_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab2_CUDA,
                        thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                        thrust::device_ptr<double> FD_Rat_CUDA,
                        thrust::device_ptr<double> D_Rat2_CUDA,
                        thrust::device_ptr<double> FSZ_Rat_CUDA,
                        thrust::device_ptr<double> SZ_Rat2_CUDA)
{
    double g2 = 2, gr, gp, gz;
    
    thrust::default_random_engine rng;

    // jump past the numbers used by the subsequences before me
    //rng.discard(N * thread_id);

    // create a mapping from random numbers to [0,1)
    thrust::uniform_real_distribution<float> u01(0,1);//printf("%f  %f %f \n", u01(rng), u01(rng), u01(rng));

    while (g2 > 1)
    {
        gr = 8*(u01(rng)-0.5);
        gp = 16*(u01(rng)-0.5);
        gz = 8*(u01(rng)-0.5);
        g2 = gr*gr/16+gp*gp/64+gz*gz/16;//cout << g2 << endl;
    }
        
    double vr = v_cyl[0]*gr;
    double vp = v_phi_max + v_cyl[1]*gp;
    double vz = v_cyl[2]*gz;
    
    //The original code defined FDisk function that simply called DiskDF5ez.
    //So I'm skipping it.
    double r_cyl = sqrt(pos[0]*pos[0]+pos[1]*pos[1]), z = pos[2];
    
    df = DiskDF5ezCUDA(vr, vp, vz, r_cyl, z, j, Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA,
                       Omega2_CUDA, A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA,
                       R_Tab2_CUDA, R_Tab2_Zero_CUDA, FD_Rat_CUDA,
                       D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA)*
                       v_cyl[0]*v_cyl[1]*v_cyl[2];
    
    double u, v, u_prime, v_prime;
    
    if (r_cyl != 0)
    {
        u = -vr*pos[0]/r_cyl + vp*pos[1]/r_cyl;
        v = -vr*pos[1]/r_cyl - vp*pos[0]/r_cyl;
    }
    else
    {
        u = (vp-vr)*0.5;
        v = -(vp+vr)*0.5;
    }
    
    RotateCoordinatesCUDA(u, v, vz, u_prime, v_prime, w_los);
}

__host__ __device__     
void DiskDispersionCUDA(double &xp, double &yp, double &wbar_disk, 
                        double &disk_dispersion, int &j,
                        thrust::device_ptr<double> Radius_CUDA,
                        thrust::device_ptr<double> Rad_Spline_CUDA,
                        thrust::device_ptr<double> A_Pot_CUDA,
                        thrust::device_ptr<double> F_R_CUDA,
                        thrust::device_ptr<double> Omega_CUDA,
                        thrust::device_ptr<double> Omega2_CUDA,
                        thrust::device_ptr<double> A_K_CUDA,
                        thrust::device_ptr<double> A_K2_CUDA,
                        thrust::device_ptr<double> Am_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab_CUDA,
                        thrust::device_ptr<double> R_Tab2_CUDA,
                        thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                        thrust::device_ptr<double> FD_Rat_CUDA,
                        thrust::device_ptr<double> D_Rat2_CUDA,
                        thrust::device_ptr<double> FSZ_Rat_CUDA,
                        thrust::device_ptr<double> SZ_Rat2_CUDA)
{
    double x, y, z, r_cyl;
    double sum = 0;
    double zdc = *(Z_Disk+j)/cos_inclination_d;
    //double rp = xp*xp+yp*yp, zdc = *(Z_Disk+j)/cos_inclination_d;
    //double dt = 2*zdc/n_los_d;
    
    double zoffset = xp*sin_inclination_d/cos_inclination_d;
    
    double zp = 0;
    
    wbar_disk = 0;
    disk_dispersion = 0;   
    
    for (int i = 1; i < n_los_d; ++i)
    {
        double t = 2.0*i/n_los_d-1;
        zp = zoffset+0.5*zdc*log((1+t)/(1-t));
        
        RotateCoordinatesCUDA(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        
        double disk_density = DiskDensfICUDA(r_cyl, z, j, Radius_CUDA, A_Pot_CUDA);
        double fac = cosh(zp/zdc)*cosh(zp/zdc);
        
        double w0 = 0, w1 = 0, w2 = 0, w_los, dfn;
        double v_cyl[3], pos[3] = {x,y,z}, f_max, v_circ, v_phi_max;
        
        PreDiskVelocitiesCUDA(pos, v_phi_max, v_cyl, f_max, v_circ, j, 
                              Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, F_R_CUDA,
                              Omega_CUDA, Omega2_CUDA, A_K_CUDA, A_K2_CUDA, 
                              Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA, 
                              R_Tab2_Zero_CUDA, FD_Rat_CUDA,
                              D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
        
        for (int k = 0; k < n_vel_d; ++k)
        {
            DiskVelocitiesCUDA(pos, v_phi_max, v_cyl, f_max, v_circ, w_los, dfn,
                               j, Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, 
                               Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
                               A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA, 
                               R_Tab2_Zero_CUDA, FD_Rat_CUDA,
                               D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
            w0 += dfn;
            w1 += dfn*w_los;
            w2 += dfn*w_los*w_los;
        }
        
        w1 /= w0;
        w2 /= w0;
        
        sum += disk_density*fac;
        wbar_disk += w1*disk_density*fac;
        disk_dispersion += w2*disk_density*fac;
    }
    
//     thrust::device_vector<double> Index(n_los_d);
//     thrust::sequence(Index.begin(), Index.end());
//     thrust::device_vector<double> Sum(n_los_d, 0);
//     thrust::device_vector<double> Wbar_Disk(n_los_d, 0);
//     thrust::device_vector<double> Disk_Disp(n_los_d, 0);
//     
//     thrust::for_each(make_zip_iterator(make_tuple(Index.begin(), Sum.begin(), 
//                                                   Wbar_Disk.begin(), Disk_Disp.begin())),
//                      make_zip_iterator(make_tuple(Index.end(), Sum.end(), 
//                                                   Wbar_Disk.end(), Disk_Disp.end())),
//                      StarDispFunctor(xp, yp, j, Radius_CUDA, Rad_Spline_CUDA,
//                                      A_Pot_CUDA, F_R_CUDA, Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
//                                      A_K2_CUDA,  
//                                      Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA, R_Tab2_Zero_CUDA,  
//                                      FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA));
//     double sum;
//     sum = thrust::reduce(Sum.begin(), Sum.end(), 0, thrust::plus<double>());
//     wbar_disk = thrust::reduce(Wbar_Disk.begin(), Wbar_Disk.end(), 0, thrust::plus<double>());
//     disk_dispersion = thrust::reduce(Disk_Disp.begin(), Disk_Disp.end(), 0, thrust::plus<double>());
    
    wbar_disk /= sum;
    disk_dispersion /= sum;
}

__host__ __device__     
void BulgeVelocitiesCUDA(double &psi, double &vmag, double &dvmag, double &x, 
                         double &y, double &z, double &w_los, double &dfn,
                         thrust::device_ptr<double> DF_Sersic_CUDA)
{
    thrust::default_random_engine rng;

    // jump past the numbers used by the subsequences before me
    //rng.discard(N * thread_id);

    // create a mapping from random numbers to [0,1)
    thrust::uniform_real_distribution<float> u01(0,1);

    double vr, vp;
    double r_cyl = sqrt(x*x+y*y);
    double cos_theta = -1+2*(u01(rng));
    double r_phi = 2*PI*(u01(rng));
    
    double u = vmag*sqrt(1-cos_theta*cos_theta)*cos(r_phi);
    double v = vmag*sqrt(1-cos_theta*cos_theta)*sin(r_phi);
    double w = vmag*cos_theta;//cout<<u<<" "<<v<<" "<<w<<" "<<endl;
    
    double E = psi - 0.5*vmag*vmag;
    dfn = vmag*vmag*vmag*dvmag*BulgeDFCUDA(E, DF_Sersic_CUDA);
    
    if (r_cyl == 0)
    {
        vr = u;
        vp = v;
    }
    else
    {
        vr = (u*x+v*y)/r_cyl;
        vp = (-u*y+v*x)/r_cyl;
    }
    
    if (u01(rng) < bulge_stream_d)
    {
        vp = -fabs(vp);
    }
    else
    {
        vp = fabs(vp);
    }
    
    if (r_cyl == 0)
    {
        u = vr;
        v = vp;
    }
    else
    {
        u = (vr*x-vp*y)/r_cyl;
        v = (vr*y+vp*x)/r_cyl;
    }
    
    double uprime, vprime;
    RotateCoordinatesCUDA(u, v, w, uprime, vprime, w_los);
}

__host__ __device__     
void BulgeDispersionCUDA(double &xp, double &yp, double &wbar_bulge, 
                         double &bulge_dispersion,
                         thrust::device_ptr<double> Radius_CUDA,
                         thrust::device_ptr<double> A_Pot_CUDA, 
                         thrust::device_ptr<double> Dens_Psi_Bulge_CUDA,
                         thrust::device_ptr<double> DF_Sersic_CUDA)
{
    double x, y, z;
    double sum = 0;
    double rp = sqrt(xp*xp+yp*yp);
    double rlogmax = log(min(a_bulge_d*pow(10.0, n_sersic_d), r_max_d)), rlogmin;
  
    if (rp > 0)
    {
        rlogmin = log(rp/100);
    }
    else 
    {
        rlogmin = log(0.1*dr_d);
    }
    
    double drlog = (rlogmax-rlogmin)/n_los_d;
    double deltar = exp(rlogmin);
    
    //First integrate over velocity with zp = 0, then zp = rlogmin
    double zp = 0;
    
    RotateCoordinatesBackCUDA(xp, yp, zp, x, y, z);
    double r_cyl = sqrt(x*x+y*y);
    double psi = PotCUDA(r_cyl, z, Radius_CUDA, A_Pot_CUDA);
    double bulge_density = BulgeDensPsiCUDA(psi, Dens_Psi_Bulge_CUDA);
    
    double w0 = 0, w1 = 0, w2 = 0, w_los, dfn;
    double v_cyl[3], v_phi_max;//, f_max, v_circ, v_phi_max;
    
    double v_max2 = 2*psi;
    double v_max = sqrt(v_max2), v_min = v_bulge_d/100;
    double vlogmax = log(v_max), vlogmin = log(v_min);
    
    double dvlog = vlogmax - vlogmin;
    
    for (int j = 0; j < n_vel_d; ++j)
    {
        double vmag = v_min*exp(dvlog*(j+1)/n_vel_d);
    
        BulgeVelocitiesCUDA(psi, vmag, dvlog, x, y, z, w_los, dfn, DF_Sersic_CUDA);
        w0 += dfn;
        w1 += dfn*w_los;
        w2 += dfn*w_los*w_los;
    }
    
    if (w0 > 0)
    {
        w1 /= w0;
        w2 /= w0;
        sum += deltar*bulge_density;
        wbar_bulge += deltar*bulge_density*w1;
        bulge_dispersion = deltar*bulge_density*w2;
    }
    
    zp = exp(rlogmin);
    
    RotateCoordinatesBackCUDA(xp, yp, zp, x, y, z);
    r_cyl = sqrt(x*x+y*y);
    psi = PotCUDA(r_cyl, z, Radius_CUDA, A_Pot_CUDA);
    bulge_density = BulgeDensPsiCUDA(psi, Dens_Psi_Bulge_CUDA);
    
    w0 = 0;
    w1 = 0;
    w2 = 0;
 
    v_max2 = 2*psi;
    v_max = sqrt(v_max2), v_min = v_bulge_d/100;
    vlogmax = log(v_max), vlogmin = log(v_min);
    
    dvlog = vlogmax - vlogmin;
    
    for (int j = 0; j < n_vel_d; ++j)
    {
        double vmag = v_min*exp(dvlog*(j+1)/n_vel_d);
   
        BulgeVelocitiesCUDA(psi, vmag, dvlog, x, y, z, w_los, dfn, DF_Sersic_CUDA);
        w0 += dfn;
        w1 += dfn*w_los;
        w2 += dfn*w_los*w_los;
   }
    
    if (w0 > 0)
    {
        w1 /= w0;
        w2 /= w0;
        sum += 0.5*(deltar+drlog*zp)*bulge_density;
        wbar_bulge += 0.5*(deltar+drlog*zp)*bulge_density*w1;
        bulge_dispersion = 0.5*(deltar+drlog*zp)*bulge_density*w2;
    }
    
    //Now do the rest of the points along the LOS
    for (int i = 2; i < n_los_d; ++i)
    {
        double zlog = rlogmin + i*drlog;
        zp = exp(zlog);
        RotateCoordinatesBackCUDA(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        psi = PotCUDA(r_cyl, z, Radius_CUDA, A_Pot_CUDA);
        bulge_density = BulgeDensPsiCUDA(psi, Dens_Psi_Bulge_CUDA);
    
        w0 = 0;
        w1 = 0;
        w2 = 0;
 
        v_max2 = 2*psi;
        v_max = sqrt(v_max2), v_min = v_bulge_d/100;
        vlogmax = log(v_max), vlogmin = log(v_min);
    
        dvlog = vlogmax - vlogmin;
    
        for (int j = 0; j < n_vel_d; ++j)
        {
            double vmag = v_min*exp(dvlog*(j+1)/n_vel_d);
        
            BulgeVelocitiesCUDA(psi, vmag, dvlog, x, y, z, w_los, dfn, DF_Sersic_CUDA);
            w0 += dfn;
            w1 += dfn*w_los;
            w2 += dfn*w_los*w_los;
        }
    
        if (w0 > 0)
        {
            w1 /= w0;
            w2 /= w0;
            sum += drlog*zp*bulge_density;
            wbar_bulge += drlog*zp*bulge_density*w1;
            bulge_dispersion += drlog*zp*bulge_density*w2;
        }
    }
    
    wbar_bulge /= sum;
    bulge_dispersion /= sum;
}

__host__ __device__     
void FindMaxCUDA(double &r, double &z, double &vsigp, double &vpmax, 
                 double &fmax, int &j,
                 thrust::device_ptr<double> Radius_CUDA,
                 thrust::device_ptr<double> Rad_Spline_CUDA,
                 thrust::device_ptr<double> A_Pot_CUDA,
                 thrust::device_ptr<double> Omega_CUDA,
                 thrust::device_ptr<double> Omega2_CUDA,
                 thrust::device_ptr<double> A_K_CUDA,
                 thrust::device_ptr<double> A_K2_CUDA,
                 thrust::device_ptr<double> Am_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab2_CUDA,
                 thrust::device_ptr<double> R_Tab2_Zero_CUDA,
                 thrust::device_ptr<double> FD_Rat_CUDA,
                 thrust::device_ptr<double> D_Rat2_CUDA,
                 thrust::device_ptr<double> FSZ_Rat_CUDA,
                 thrust::device_ptr<double> SZ_Rat2_CUDA)
{
    double dv = 0.1*vsigp, vpm = vpmax, vpmold, zero=0;
    double v0 = vpm-dv, v1 = vpm+dv;
    double f0 = DiskDF5ezCUDA(zero,v0,zero,r,z,j,
                              Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA, Omega2_CUDA, 
                              A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, 
                              R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                              FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
    double fmid = DiskDF5ezCUDA(zero,vpm,zero,r,z,j,
                                Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA, Omega2_CUDA, 
                                A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, 
                                R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                                FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
    double f1 = DiskDF5ezCUDA(zero,v1,zero,r,z,j,
                              Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA, Omega2_CUDA, 
                              A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, 
                              R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                              FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);

    if (fmid>=f0 && fmid>f1)
    {
        fmax = fmid;
    }
    else
    {
        if (f0 > f1)
        {
            double ftmp = f0;
            f0 = f1;
            f1 = ftmp;
            v1 = v0;
            dv = -dv;
        }
        
        vpm = v1;
        
        int flag = 1;
        
        while(flag > 0)
        {
            dv *= 2;
            vpmold = vpm;
            vpm += dv;
            f0 = f1;
            f1 = DiskDF5ezCUDA(zero,vpm,zero,r,z,j,
                               Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA, Omega2_CUDA, 
                               A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, 
                               R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                               FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
            flag = (f1>f0);
        }
        
        vpmax = vpmold;
        fmax = f0;
    }
}            

