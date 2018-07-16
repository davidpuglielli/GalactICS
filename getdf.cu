struct GetNFWDFFunctor
{
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> H_Pot_CUDA;
    thrust::device_ptr<double> B_Pot_CUDA;
    thrust::device_ptr<double> D_Pot_CUDA;
    thrust::device_ptr<double> H_FR_CUDA;
    thrust::device_ptr<double> B_FR_CUDA;
    thrust::device_ptr<double> D_FR_CUDA;
    thrust::device_ptr<double> D_Dens_CUDA;
    thrust::device_ptr<double> Halo_AC_Radius_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_D_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA;
    
    GetNFWDFFunctor(thrust::device_ptr<double> _Radius_CUDA,
                    thrust::device_ptr<double> _H_Pot_CUDA,
                    thrust::device_ptr<double> _B_Pot_CUDA,
                    thrust::device_ptr<double> _D_Pot_CUDA,
                    thrust::device_ptr<double> _H_FR_CUDA,
                    thrust::device_ptr<double> _B_FR_CUDA,
                    thrust::device_ptr<double> _D_FR_CUDA,
                    thrust::device_ptr<double> _D_Dens_CUDA,
                    thrust::device_ptr<double> _Halo_AC_Radius_CUDA,
                    thrust::device_ptr<double> _Halo_AC_Dens_CUDA,
                    thrust::device_ptr<double> _Halo_AC_Dens_D_CUDA,
                    thrust::device_ptr<double> _Halo_AC_Dens_DD_CUDA):
                    Radius_CUDA(_Radius_CUDA),
                    H_Pot_CUDA(_H_Pot_CUDA),
                    B_Pot_CUDA(_B_Pot_CUDA),
                    D_Pot_CUDA(_D_Pot_CUDA),
                    H_FR_CUDA(_H_FR_CUDA),
                    B_FR_CUDA(_B_FR_CUDA),
                    D_FR_CUDA(_D_FR_CUDA),
                    D_Dens_CUDA(_D_Dens_CUDA),
                    Halo_AC_Radius_CUDA(_Halo_AC_Radius_CUDA),
                    Halo_AC_Dens_CUDA(_Halo_AC_Dens_CUDA),
                    Halo_AC_Dens_D_CUDA(_Halo_AC_Dens_D_CUDA),
                    Halo_AC_Dens_DD_CUDA(_Halo_AC_Dens_DD_CUDA){}
    
    __host__ __device__
    float operator()(double &energy) const
    {
        //printf("in functor\n");
        
        double t, t_max = sqrt(energy-psi_crit_d), dt = t_max/(n_int_d - 1);
        double psi = energy;
        double sum = 0, r_upper, r_lower;

        FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
                         B_Pot_CUDA, D_Pot_CUDA);
        double tolerance = tolerance_factor_d;
        double r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, 
                                         Radius_CUDA, H_Pot_CUDA, 
                                         B_Pot_CUDA, D_Pot_CUDA);

        double d2rho_dpsi2 = GetNFWd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
                                                  B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA,
                                                  Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA,
                                                  Halo_AC_Dens_D_CUDA, Halo_AC_Dens_DD_CUDA);
        sum += dt*d2rho_dpsi2;

        for (int i = 0; i < n_int_d-1; ++i)
        {
            t = (i + 1)*dt;
            psi = energy -t*t;
            FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
            tolerance = tolerance_factor_d;
            r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
            d2rho_dpsi2 = GetNFWd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
                                               B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA,
                                               Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA,
                                               Halo_AC_Dens_D_CUDA, Halo_AC_Dens_DD_CUDA);

            sum += 2*dt*d2rho_dpsi2;
        }

        //printf("in functor2\n");
        return sum*oneoversqrt8pi2;
    }
    
//     __host__ __device__
//     float operator()(double &energy, int &index) const
//     {
//         //printf("in functor %d\n", index);
//         if (index==n_int_d) return 0;
//         int i = index;
//         
//         double t_max = sqrt(energy-psi_crit_d), dt = t_max/(n_int_d - 1);
//         double t, psi;
//         double sum = 0, r_upper, r_lower;
//         
//         t = i*dt;
//         psi = energy-t*t;
// 
//         FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
//                          B_Pot_CUDA, D_Pot_CUDA);
//         double tolerance = tolerance_factor_d;
//         double r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, 
//                                          Radius_CUDA, H_Pot_CUDA, 
//                                          B_Pot_CUDA, D_Pot_CUDA);
// 
//         double d2rho_dpsi2 = GetNFWd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
//                                                   B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA);
//         
//         if (index>0)
//         {
//             return 2*dt*d2rho_dpsi2;
//         }
//         else
//         {
//             return dt*d2rho_dpsi2;
//         }
//     }
};

struct GetSersicDFFunctor
{
    thrust::device_ptr<double> Radius_CUDA;
    thrust::device_ptr<double> H_Pot_CUDA;
    thrust::device_ptr<double> B_Pot_CUDA;
    thrust::device_ptr<double> D_Pot_CUDA;
    thrust::device_ptr<double> H_FR_CUDA;
    thrust::device_ptr<double> B_FR_CUDA;
    thrust::device_ptr<double> D_FR_CUDA;
    thrust::device_ptr<double> D_Dens_CUDA;
    thrust::device_ptr<double> Halo_AC_Radius_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_D_CUDA;
    thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA;
    
    GetSersicDFFunctor(thrust::device_ptr<double> _Radius_CUDA,
                       thrust::device_ptr<double> _H_Pot_CUDA,
                       thrust::device_ptr<double> _B_Pot_CUDA,
                       thrust::device_ptr<double> _D_Pot_CUDA,
                       thrust::device_ptr<double> _H_FR_CUDA,
                       thrust::device_ptr<double> _B_FR_CUDA,
                       thrust::device_ptr<double> _D_FR_CUDA,
                       thrust::device_ptr<double> _D_Dens_CUDA,
                       thrust::device_ptr<double> _Halo_AC_Radius_CUDA,
                       thrust::device_ptr<double> _Halo_AC_Dens_CUDA,
                       thrust::device_ptr<double> _Halo_AC_Dens_D_CUDA,
                       thrust::device_ptr<double> _Halo_AC_Dens_DD_CUDA):
                       Radius_CUDA(_Radius_CUDA),
                       H_Pot_CUDA(_H_Pot_CUDA),
                       B_Pot_CUDA(_B_Pot_CUDA),
                       D_Pot_CUDA(_D_Pot_CUDA),
                       H_FR_CUDA(_H_FR_CUDA),
                       B_FR_CUDA(_B_FR_CUDA),
                       D_FR_CUDA(_D_FR_CUDA),
                       D_Dens_CUDA(_D_Dens_CUDA),
                       Halo_AC_Radius_CUDA(_Halo_AC_Radius_CUDA),
                       Halo_AC_Dens_CUDA(_Halo_AC_Dens_CUDA),
                       Halo_AC_Dens_D_CUDA(_Halo_AC_Dens_D_CUDA),
                       Halo_AC_Dens_DD_CUDA(_Halo_AC_Dens_DD_CUDA){}
    
    __host__ __device__  
    float operator()(double &energy) const
    {
        //printf("in functor\n");
        
        double t, t_max = sqrt(energy-psi_crit_d), dt = t_max/(n_int_d - 1);
        double psi = energy;
        double sum = 0, r_upper, r_lower;

        //printf("in functor60\n");
        FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);

        double tolerance = tolerance_factor_d;
        double r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, 
                                         Radius_CUDA, H_Pot_CUDA, B_Pot_CUDA, D_Pot_CUDA);
        //printf("in functor %.12f %.12f %.12f %.12f\n", psi, r_psi, r_lower, r_upper);

        //printf("in functor6\n");
        double d2rho_dpsi2 = GetSersicd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
                                                     B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA,
                                                  Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA,
                                                  Halo_AC_Dens_D_CUDA, Halo_AC_Dens_DD_CUDA);
        sum += dt*d2rho_dpsi2;
        //printf("in functor %.12f %.12f %.12f %.12f\n", r_psi, psi, sum, d2rho_dpsi2);

        //printf("in functor1\n");
        
        for (int i = 0; i < n_int_d-1; ++i)
        {
            t = (i + 1)*dt;
            psi = energy -t*t;
            FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
            tolerance = tolerance_factor_d;
            r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
            d2rho_dpsi2 = GetSersicd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
                                                  B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA,
                                                  Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA,
                                                  Halo_AC_Dens_D_CUDA, Halo_AC_Dens_DD_CUDA);

            sum += 2*dt*d2rho_dpsi2;
        }

        //printf("in functor2\n");
        
        return sum*oneoversqrt8pi2;
    }
    
//     __host__ __device__
//     float operator()(double &energy, int &index) const
//     {
//         //printf("in functor\n");
//         if (index==n_int_d-1) return 0;
//         int i = index;
//         
//         double t_max = sqrt(energy-psi_crit_d), dt = t_max/(n_int_d - 1);
//         double t, psi;
//         double sum = 0, r_upper, r_lower;
//         
//         t = i*dt;
//         psi = energy-t*t;
// 
//         FindBracketsCUDA(psi, r_lower, r_upper, Radius_CUDA, H_Pot_CUDA, 
//                          B_Pot_CUDA, D_Pot_CUDA);
//         double tolerance = tolerance_factor_d;
//         double r_psi = RootBisectionCUDA(psi, r_upper, r_lower, tolerance, 
//                                          Radius_CUDA, H_Pot_CUDA, B_Pot_CUDA, D_Pot_CUDA);
//         //printf("in functor %d %.12f %.12f %.12f %.12f\n", index, psi, r_psi, r_lower, r_upper);
// 
//         double d2rho_dpsi2 = GetSersicd2rhodpsi2CUDA(r_psi, Radius_CUDA, H_FR_CUDA, 
//                                                      B_FR_CUDA, D_FR_CUDA, D_Dens_CUDA);
//         
//         //if(psi>33)printf("in functor %d %.12f %.12f %.12f %.12f %.12f\n", 
//         //                   index, psi, r_psi, r_lower, r_upper, d2rho_dpsi2);
// 
//         if (index>0)
//         {
//             return 2*dt*d2rho_dpsi2;
//         }
//         else
//         {
//             return dt*d2rho_dpsi2;
//         }
//     }
};

void GenNFWDistFuncCUDA(void)
{
    //thrust::device_ptr<double> Table_E_CUDA = thrust::device_malloc<double>(n_psi);
    thrust::device_ptr<double> Table_E_CUDA = thrust::device_malloc<double>(n_psi);//*n_int);
    //thrust::device_ptr<int> Index_CUDA = thrust::device_malloc<int>(n_psi*n_int);
    thrust::device_ptr<double> DF_NFW_CUDA = thrust::device_malloc<double>(n_psi);//*n_int);
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> H_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> B_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> H_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> B_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_Dens_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Halo_AC_Radius_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_D_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA = thrust::device_malloc<double>(nr_ac);
    
    thrust::copy(Table_E.begin(), Table_E.end(), Table_E_CUDA);
//     for (int i = 0; i < n_psi; ++i)
//     {
//         int offset = i*n_int, offset2 = (i+1)*n_int;
//         double table_e = Table_E[i];
//         thrust::fill(Table_E_CUDA+offset, Table_E_CUDA+offset2, table_e);
//         thrust::sequence(Index_CUDA+offset, Index_CUDA+offset2);
//     }
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    thrust::copy(H_Pot.begin(), H_Pot.end(), H_Pot_CUDA);
    thrust::copy(B_Pot.begin(), B_Pot.end(), B_Pot_CUDA);
    thrust::copy(D_Pot.begin(), D_Pot.end(), D_Pot_CUDA);
    thrust::copy(H_FR.begin(), H_FR.end(), H_FR_CUDA);
    thrust::copy(B_FR.begin(), B_FR.end(), B_FR_CUDA);
    thrust::copy(D_FR.begin(), D_FR.end(), D_FR_CUDA);
    thrust::copy(D_Dens.begin(), D_Dens.end(), D_Dens_CUDA);
    thrust::copy(Halo_AC_Radius.begin(), Halo_AC_Radius.end(), Halo_AC_Radius_CUDA);
    thrust::copy(Halo_AC_Dens.begin(), Halo_AC_Dens.end(), Halo_AC_Dens_CUDA);
    thrust::copy(Halo_AC_Dens_D.begin(), Halo_AC_Dens_D.end(), Halo_AC_Dens_D_CUDA);
    thrust::copy(Halo_AC_Dens_DD.begin(), Halo_AC_Dens_DD.end(), Halo_AC_Dens_DD_CUDA);
    //thrust::copy(DF_NFW.begin(), DF_NFW.end(), DF_NFW_CUDA);
    
    //int size = n_psi*n_int;
    thrust::transform(Table_E_CUDA, Table_E_CUDA+n_psi, DF_NFW_CUDA, 
                      GetNFWDFFunctor(Radius_CUDA, H_Pot_CUDA, B_Pot_CUDA, 
                                      D_Pot_CUDA, H_FR_CUDA, B_FR_CUDA, 
                                      D_FR_CUDA, D_Dens_CUDA, Halo_AC_Radius_CUDA,
                                      Halo_AC_Dens_CUDA, Halo_AC_Dens_D_CUDA,
                                      Halo_AC_Dens_DD_CUDA));
    
//     for (int i = 0; i < n_psi; ++i)
//     {
//         int offset = i*n_int, offset2 = (i+1)*n_int;
//         double df_nfw = 0;
//         for (int j=0; j < n_int; ++j)
//         {
//             df_nfw += *(DF_NFW_CUDA+offset+j);
//         }
//         //double df_nfw = thrust::reduce(DF_NFW_CUDA+offset, DF_NFW_CUDA+offset2, 
//         //                               0, thrust::plus<double>());
//         DF_NFW[i] = df_nfw*oneoversqrt8pi2;
//         //cout << i << 
//     }
    
    thrust::copy(DF_NFW_CUDA, DF_NFW_CUDA+n_psi, DF_NFW.begin());
    
    double DF_NFW_last = -1000;
    
    ofstream DFFile;
    
    if (do_file_io)
    {
        DFFile.open("dfnfw.dat", ios::out);
    }
    
    for (int i = 0; i < n_psi; ++i)
    {
        if (DF_NFW[i] > 0)
        {
            DF_NFW[i] = log(DF_NFW[i]);
            DF_NFW_last = DF_NFW[i];
        }
        else
        {
            cout << "Warning: DF_NFW < 0 at bin " << i << ". Using previous "
                 << "value. " << DF_NFW[i] << endl;
            DF_NFW[i] = DF_NFW_last;
        }
        
        if (do_file_io)
        {
            DFFile << setprecision(12) << Table_E[i] << "  " << DF_NFW[i] << endl;
        }
    }
    
    thrust::device_free(Table_E_CUDA);
    thrust::device_free(DF_NFW_CUDA);
    thrust::device_free(Radius_CUDA);
    thrust::device_free(H_Pot_CUDA);
    thrust::device_free(B_Pot_CUDA);
    thrust::device_free(D_Pot_CUDA);
    thrust::device_free(H_FR_CUDA);
    thrust::device_free(B_FR_CUDA);
    thrust::device_free(D_FR_CUDA);
    thrust::device_free(D_Dens_CUDA);
    thrust::device_free(Halo_AC_Radius_CUDA);
    thrust::device_free(Halo_AC_Dens_CUDA);
    thrust::device_free(Halo_AC_Dens_D_CUDA);
    thrust::device_free(Halo_AC_Dens_DD_CUDA);
}

void GenSersicDistFuncCUDA(void)
{
    thrust::device_ptr<double> Table_E_CUDA = thrust::device_malloc<double>(n_psi);//*n_int);
    //thrust::device_ptr<int> Index_CUDA = thrust::device_malloc<int>(n_psi*n_int);
    thrust::device_ptr<double> DF_Sersic_CUDA = thrust::device_malloc<double>(n_psi);//*n_int);
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> H_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> B_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_Pot_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> H_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> B_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_FR_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_Dens_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Halo_AC_Radius_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_D_CUDA = thrust::device_malloc<double>(nr_ac);
    thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA = thrust::device_malloc<double>(nr_ac);
    
    thrust::copy(Table_E.begin(), Table_E.end(), Table_E_CUDA);
//     for (int i = 0; i < n_psi; ++i)
//     {
//         int offset = i*n_int, offset2 = (i+1)*n_int;
//         double table_e = Table_E[i];
//         //cout << i << " " << offset << " " << offset2 << " " << table_e << endl;
//         thrust::fill(Table_E_CUDA+offset, Table_E_CUDA+offset2, table_e);
//         thrust::sequence(Index_CUDA+offset, Index_CUDA+offset2);
//     }
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    thrust::copy(H_Pot.begin(), H_Pot.end(), H_Pot_CUDA);
    thrust::copy(B_Pot.begin(), B_Pot.end(), B_Pot_CUDA);
    thrust::copy(D_Pot.begin(), D_Pot.end(), D_Pot_CUDA);
    thrust::copy(H_FR.begin(), H_FR.end(), H_FR_CUDA);
    thrust::copy(B_FR.begin(), B_FR.end(), B_FR_CUDA);
    thrust::copy(D_FR.begin(), D_FR.end(), D_FR_CUDA);
    thrust::copy(D_Dens.begin(), D_Dens.end(), D_Dens_CUDA);
    thrust::copy(Halo_AC_Radius.begin(), Halo_AC_Radius.end(), Halo_AC_Radius_CUDA);
    thrust::copy(Halo_AC_Dens.begin(), Halo_AC_Dens.end(), Halo_AC_Dens_CUDA);
    thrust::copy(Halo_AC_Dens_D.begin(), Halo_AC_Dens_D.end(), Halo_AC_Dens_D_CUDA);
    thrust::copy(Halo_AC_Dens_DD.begin(), Halo_AC_Dens_DD.end(), Halo_AC_Dens_DD_CUDA);
    
    //int size = n_psi*n_int;
    thrust::transform(Table_E_CUDA, Table_E_CUDA+n_psi, DF_Sersic_CUDA, 
                      GetSersicDFFunctor(Radius_CUDA, H_Pot_CUDA, B_Pot_CUDA, 
                                         D_Pot_CUDA, H_FR_CUDA, B_FR_CUDA, D_FR_CUDA, 
                                         D_Dens_CUDA, Halo_AC_Radius_CUDA,
                                         Halo_AC_Dens_CUDA, Halo_AC_Dens_D_CUDA,
                                         Halo_AC_Dens_DD_CUDA));
    
//     for (int i = 0; i < n_psi; ++i)
//     {
//         int offset = i*n_int, offset2 = (i+1)*n_int;
//         double df_sersic = 0;
//         for (int j=0; j < n_int; ++j)
//         {
//             df_sersic += *(DF_Sersic_CUDA+offset+j);
//         }
//         //double df_sersic = thrust::reduce(DF_Sersic_CUDA+offset, DF_Sersic_CUDA+offset2, 
//         //                                  0, thrust::plus<double>());
//         DF_Sersic[i] = df_sersic*oneoversqrt8pi2;
//         //cout << i << " " << df_sersic*oneoversqrt8pi2 << endl;
//     }
    
    thrust::copy(DF_Sersic_CUDA, DF_Sersic_CUDA+n_psi, DF_Sersic.begin());
    
    double DF_Sersic_last = -1000;
    
    ofstream DFFile;
    
    if (do_file_io)
    {
        DFFile.open("dfsersic.dat", ios::out);
    }
    
    for (int i = 0; i < n_psi; ++i)
    {
        if (DF_Sersic[i] > 0)
        {
            DF_Sersic[i] = log(DF_Sersic[i]);
            DF_Sersic_last = DF_Sersic[i];
        }
        else
        {
            cout << "Warning: DF_Sersic < 0 at bin " << i << ". Using previous "
                 << "value. " << DF_Sersic[i] << endl;
            DF_Sersic[i] = DF_Sersic_last;
        }
        
        if (do_file_io)
        {
            DFFile << setprecision(12) << Table_E[i] << "  " << DF_Sersic[i] << endl;
        }
    }
    
    thrust::device_free(Table_E_CUDA);
    thrust::device_free(DF_Sersic_CUDA);
    thrust::device_free(Radius_CUDA);
    thrust::device_free(H_Pot_CUDA);
    thrust::device_free(B_Pot_CUDA);
    thrust::device_free(D_Pot_CUDA);
    thrust::device_free(H_FR_CUDA);
    thrust::device_free(B_FR_CUDA);
    thrust::device_free(D_FR_CUDA);
    thrust::device_free(D_Dens_CUDA);
    thrust::device_free(Halo_AC_Radius_CUDA);
    thrust::device_free(Halo_AC_Dens_CUDA);
    thrust::device_free(Halo_AC_Dens_D_CUDA);
    thrust::device_free(Halo_AC_Dens_DD_CUDA);
}

//what about the cusp?
//Wasteful - we calculate the force twice when we should only need to once. Fix!
__host__ __device__   
double GetNFWd2rhodpsi2CUDA(double &r, 
                            thrust::device_ptr<double> Radius_CUDA,
                            thrust::device_ptr<double> H_FR_CUDA,
                            thrust::device_ptr<double> B_FR_CUDA,
                            thrust::device_ptr<double> D_FR_CUDA,
                            thrust::device_ptr<double> D_Dens_CUDA,
                            thrust::device_ptr<double> Halo_AC_Radius_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_D_CUDA,
                            thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA)
{
    double density = HaloDensityCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA);
    double ddensity = HaloDensityPrimeCUDA(r, Halo_AC_Radius_CUDA, 
                                           Halo_AC_Dens_CUDA, Halo_AC_Dens_D_CUDA);
    double dddensity = HaloDensity2PrimeCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA,
                                             Halo_AC_Dens_D_CUDA, Halo_AC_Dens_DD_CUDA);
    double force = HaloForceCUDA(r, Radius_CUDA, H_FR_CUDA);
    
    
    if (bulge_flag_d == 1)
    {
        force += SersicForceCUDA(r, Radius_CUDA, B_FR_CUDA);
        density += SersicDensCUDA(r);
    }
    
    if (disk_flag_d == 1)
    {
        force += DiskForceCUDA(r, Radius_CUDA, D_FR_CUDA);
        density += DiskDensityCUDA(r, Radius_CUDA, D_Dens_CUDA);
    }
    
    double t1 = fourpi*density*ddensity/force;
    double t2 = 2*ddensity/r;

    return (t1 + t2 + dddensity)/(force*force);
}

__host__ __device__   
double GetSersicd2rhodpsi2CUDA(double &r, 
                               thrust::device_ptr<double> Radius_CUDA,
                               thrust::device_ptr<double> H_FR_CUDA,
                               thrust::device_ptr<double> B_FR_CUDA,
                               thrust::device_ptr<double> D_FR_CUDA,
                               thrust::device_ptr<double> D_Dens_CUDA,
                               thrust::device_ptr<double> Halo_AC_Radius_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_D_CUDA,
                               thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA)
{
    double density = SersicDensCUDA(r);
    double ddensity = SersicDensPrimeCUDA(r);
    double dddensity = SersicDens2PrimeCUDA(r);
    double force = SersicForceCUDA(r, Radius_CUDA, B_FR_CUDA);
    //if(r<0.01)printf("drho %.12f %.12f %.12f %.12f %.12f\n", r, density, ddensity, dddensity, force);
    
    if (halo_flag_d == 1)
    {
        force += HaloForceCUDA(r, Radius_CUDA, H_FR_CUDA);
        density += HaloDensityCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA);
    }
    //if(r<0.01)printf("drho %.12f %.12f %.12f %.12f %.12f\n", r, density, ddensity, dddensity, force);
    
    if (disk_flag_d == 1)
    {
        force += DiskForceCUDA(r, Radius_CUDA, D_FR_CUDA);
        density += DiskDensityCUDA(r, Radius_CUDA, D_Dens_CUDA);
    }
    
    //if(r<0.01)printf("drho %.12f %.12f %.12f %.12f %.12f\n", r, density, ddensity, dddensity, force);
    double t1 = fourpi*density*ddensity/force;
    double t2 = 2*ddensity/r;

    return (t1 + t2 + dddensity)/(force*force);
}

__host__ __device__   
void FindBracketsCUDA(double &psi, double &r_lower, double &r_upper, 
                      thrust::device_ptr<double> Radius_CUDA,
                      thrust::device_ptr<double> H_Pot_CUDA,
                      thrust::device_ptr<double> B_Pot_CUDA,
                      thrust::device_ptr<double> D_Pot_CUDA)
{
     r_lower = r_upper = 1;
     
     double psi_lower = psi_lower_initial_d;//GetTotalPsi(r_lower);
     double psi_upper = psi_upper_initial_d;//GetTotalPsi(r_upper);
     
     while(psi_lower <= psi)
     {
         r_lower*=0.1;
         psi_lower=GetTotalPsiCUDA(r_lower, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
     }
     
     while (psi_upper >= psi)
     {
         r_upper=r_upper*10;
         psi_upper=GetTotalPsiCUDA(r_upper, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA);
     }
}

__host__ __device__   
double RootBisectionCUDA(double &psi, double &r_upper, double &r_lower, 
                         double tolerance_d, thrust::device_ptr<double> Radius_CUDA,
                         thrust::device_ptr<double> H_Pot_CUDA,
                         thrust::device_ptr<double> B_Pot_CUDA,
                         thrust::device_ptr<double> D_Pot_CUDA)
{
    double psi_upper = GetTotalPsiCUDA(r_upper, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA)-psi;
    double psi_lower = GetTotalPsiCUDA(r_lower, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA)-psi;
    double r_midpoint, psi_midpoint, r_gap;
    
    int iters = 0;
    
    if (psi_upper*psi_lower > 0)
    {
        printf("Bisection endpoints do not bracket root. Exiting...\n %f %f %f %f\n",
               psi_upper, psi_lower, r_lower, r_upper);
    }
    
    do
    {
        ++iters;
        
        r_midpoint = (r_upper+r_lower)/2;
        
        psi_midpoint = GetTotalPsiCUDA(r_midpoint, Radius_CUDA, H_Pot_CUDA, 
                                       B_Pot_CUDA, D_Pot_CUDA)-psi;
        
        if (psi_lower*psi_midpoint < 0)
        {
            r_upper = r_midpoint;
            psi_upper = psi_midpoint;
        }
        else if (psi_upper*psi_midpoint < 0)
        {
            r_lower = r_midpoint;
            psi_lower = psi_midpoint;
        }
        else if (psi_midpoint==0)
        {
            return r_midpoint;
        }
        else
        {
            printf("Root bisection failed! Exiting...\n");
            //exit(1);
        }
        
        r_gap = r_upper-r_lower;
    }
    while (iters < max_iter_d && r_gap > tolerance_d);
    
    return r_midpoint;
}

__host__ __device__   
double GetTotalPsiCUDA(double &r, thrust::device_ptr<double> Radius_CUDA,
                       thrust::device_ptr<double> H_Pot_CUDA,
                       thrust::device_ptr<double> B_Pot_CUDA,
                       thrust::device_ptr<double> D_Pot_CUDA)
{
    double gettotalpsi = 0;
    
    //double log_r = log10(r);
    //int ihi = ceil((log_r-log_dr)/delta_logr+1);
    int ihi = ceil(r/dr_d);
    
    //if (ihi < 0) ihi = 0;
    //else if (ihi > nr-1) ihi = nr-1;
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
    double tm = 1 - t;
    
    if (disk_flag_d == 1)
    {
        gettotalpsi += t*D_Pot_CUDA[ihi] + tm*D_Pot_CUDA[ihi-1];
        //double val1 = *(D_Pot_CUDA+ihi);
        //double val2 = D_Pot_CUDA[ihi];
        //if(r<1)printf("psid %f %f %f %f\n", r, gettotalpsi, val1, val2);
    }
    //if (gasdisk_flag_d == 1)
    //{
    //    gettotalpsi += t*G_Pot_CUDA[ihi] + tm*G_Pot_CUDA[ihi-1];
    //}
    if (bulge_flag_d == 1)
    {
        if (sersic_flag_d)
        {
            gettotalpsi += SersicPotentialCUDA(r);
        }
        else
        {
            gettotalpsi += t*B_Pot_CUDA[ihi] + tm*B_Pot_CUDA[ihi-1];
            //printf("psib %f\n", gettotalpsi);
        }
    }
    if (halo_flag_d == 1)
    {
        gettotalpsi += t*H_Pot_CUDA[ihi] + tm*H_Pot_CUDA[ihi-1];
        //printf("psih %f\n", gettotalpsi);
    }
    
    return gettotalpsi;
}

__host__ __device__   
double DiskForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                     thrust::device_ptr<double> D_FR_CUDA)
{
    //double log_r = log10(r);
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
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_d-1;
    }
    
    double r1 = Radius_CUDA[ihi-1];
    double r2 = Radius_CUDA[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*D_FR_CUDA[ihi] + tm1*D_FR_CUDA[ihi-1];
}

__host__ __device__   
double DiskDensityCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                       thrust::device_ptr<double> D_Dens_CUDA)
{
    //double log_r = log10(r);
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
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_d-1;
    }
    
    double r1 = Radius_CUDA[ihi-1];
    double r2 = Radius_CUDA[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*D_Dens_CUDA[ihi] + tm1*D_Dens_CUDA[ihi-1];
}

