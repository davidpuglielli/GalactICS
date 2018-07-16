// struct SurfaceBrightnessFunctor
// {
//     const double ;
//     const double error_factor;
//     const double error_factor;
//     thrust::device_ptr<double> Radius_CUDA;
//     thrust::device_ptr<double> A_Pot_CUDA;
//     thrust::device_ptr<double> F_R_CUDA;
//     thrust::device_ptr<double> DF_Sersic_CUDA;
//     thrust::device_ptr<double> Dens_Psi_Bulge_CUDA;
// 
//     SurfaceBrightnessFunctor(double _radius, double _d_phi
//                              thrust::device_ptr<double> _Radius_CUDA,
//                              thrust::device_ptr<double> _A_Pot_CUDA,
//                              thrust::device_ptr<double> _F_R_CUDA,
//                              thrust::device_ptr<double> _DF_Sersic_CUDA,
//                              thrust::device_ptr<double> _Dens_Psi_Bulge_CUDA): 
//                              error_factor(_error_factor),
//                              Radius_CUDA(_Radius_CUDA),
//                              A_Pot_CUDA(_A_Pot_CUDA),
//                              F_R_CUDA(_F_R_CUDA),
//                              DF_Sersic_CUDA(_DF_Sersic_CUDA),
//                              Dens_Psi_Bulge_CUDA(_Dens_Psi_Bulge_CUDA){}
//     //StarVelocityFunctor(){}
// 
//     //SurfaceBrightnessFunctor(double _error_factor): error_factor(_error_factor){}
//     //SurfaceBrightnessFunctor(){}
//     
//     template <typename Tuple>
//     __host__ __device__ 
//     float operator()(const int& j) const
//     {
//             double yp = radius*cos(j*d_phi);
//             double xp = radius*sin(j*d_phi)*cos(point_inclination);
//             //cout << "step " << i << " " << j << " " << xp << " " << yp << endl;
//             
//             for (int k = 0; k < disk_d; ++k)
//             {
//                 double disk_surf_den = DiskSurfaceDensityCUDA(xp, yp, k, Radius_CUDA, A_Pot_CUDA);
//                 integral_disk[k] += disk_surf_den;
//             }
//             
//             double bulge_surf_den = BulgeSurfaceDensityCUDA(xp, yp, Radius_CUDA,
//                                                             A_Pot_CUDA, Dens_Psi_Bulge_CUDA);
//             
//             if (j == 0 || j == 99)
//             {
//                 bulge_surf_den *= 0.5;
//             }
//             
//             integral_bulge += bulge_surf_den;
//         }
//         
//         //printf("in functor3\n");
//         double bulge_light = 0, disk_light_tot = 0;
//         
//         for (int j = 0; j < disk_d; ++j)
//         {
//             integral_disk[j] /= 100;
//             disk_light[j] += integral_disk[j]/(*(ML_Disk+j))*m_scale/kpc_to_arcsec/kpc_to_arcsec;
//             disk_light_tot += disk_light[j];
//         }
//         
//         integral_bulge /= 100;
//         bulge_light = integral_bulge/ml_bulge_d*m_scale/kpc_to_arcsec/kpc_to_arcsec;
//         
//         //the 1e-6 is for conv to pc in log arg
//         double disk_magnitude_tot = mag_sun-5-2.5*log10(disk_light_tot/distance_d/distance_d*1e-6);
//         double bulge_magnitude = mag_sun-5-2.5*log10(bulge_light/distance_d/distance_d*1e-6);
//         double total_magnitude = mag_sun-5-2.5*log10((disk_light_tot+bulge_light)/distance_d/distance_d*1e-6);
//         double *disk_magnitude = new double[disk_d];
//         
//         for (int j = 0; j < disk_d; ++j)
//         {
//             disk_magnitude[j] = mag_sun-5-2.5*log10(disk_light[j]/distance_d/distance_d*1e-6);
//         }
//         
//         //printf("in functor4 %f %f %f\n", total_magnitude, data, error);
//         double chi_square = (total_magnitude-data)*(total_magnitude-data)/error/error + 2*log(error); 
//         
//         //printf("in functor5\n");
//         return chi_square;       
//     }
// };

struct StarDispFunctor
{
    const double xp;
    const double yp;
    const int j;
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

    StarDispersionFunctor(double _xp, double _yp,int _j,
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
                          xp(_xp), yp(_yp), j(_j),
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
    __host__ __device__ 
    float operator()(const Tuple data_point) const
    {
        double i = thrust::get<0>(data_point);
        
        double rp = xp*xp+yp*yp, zdc = *(Z_Disk+j)/cos_inclination_d;
        double zoffset = xp*sin_inclination_d/cos_inclination_d;
        double t = 2.0*i/n_los_d-1;
        double zp = zoffset+0.5*zdc*log((1+t)/(1-t));
        
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
        
        //Return sum, wbar and disp to transform where they will be added
        thrust::get<1>(data_point) = disk_density*fac;
        thrust::get<2>(data_point) = w1*disk_density*fac;
        thrust::get<3>(data_point) = w2*disk_density*fac;
        
    }
};
