//This file includes host functions that copy data to the device and
//call the device code through thrust::transform

#include "galaxy_cuda.h"
//#include "cuda_layer2functors.cu"
#include "cuda_functors.cu"
#include "cuda_functions.cu"
#include "getdf.cu"
#include "haloprofiles.cu"
#include "sersicprofiles.cu"
#include "diskprofiles.cu"
#include "gasdiskprofiles.cu"
#include "getp.cu"
#include "getchisquarecuda.cu"

void CopyGlobalsToDevice1(void)
{
    cout << "copy" << endl;
    
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
        cudaMemcpyToSymbol(Gamma, &G.Gamma[0], sizeof(Gamma));
        cudaMemcpyToSymbol(GasDisk_Const_d, &GasDisk_Const[0], sizeof(GasDisk_Const));
        cudaMemcpyToSymbol(Rho_GasDisk_Const_d, &Rho_GasDisk_Const[0], sizeof(Rho_GasDisk_Const));
    }
    
    cudaMemcpyToSymbol(halo_const_d, &halo_const, sizeof(double));
    cudaMemcpyToSymbol(rho_0_d, &G.rho_0, sizeof(double));
    
    cudaMemcpyToSymbol(halo_flag_d, &halo_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_flag_d, &disk_flag, sizeof(int));
    cudaMemcpyToSymbol(gasdisk_flag_d, &gasdisk_flag, sizeof(int));
    cudaMemcpyToSymbol(bulge_flag_d, &bulge_flag, sizeof(int));
    cudaMemcpyToSymbol(n_psi_d, &n_psi, sizeof(int));
    cudaMemcpyToSymbol(nr_d, &nr, sizeof(int));
    cudaMemcpyToSymbol(smbh_flag_d, &smbh_flag, sizeof(int));
    cudaMemcpyToSymbol(nbody_flag_d, &nbody_flag, sizeof(int));
    cudaMemcpyToSymbol(do_file_io_d, &do_file_io, sizeof(int));
    cudaMemcpyToSymbol(chisq_flag_d, &chisq_flag, sizeof(int));
    cudaMemcpyToSymbol(sersic_flag_d, &sersic_flag, sizeof(int));
    cudaMemcpyToSymbol(disk_d, &disk, sizeof(int));
    cudaMemcpyToSymbol(gas_disk_d, &gas_disk, sizeof(int));
    cudaMemcpyToSymbol(nr_ac_d, &nr_ac, sizeof(int));
    cudaMemcpyToSymbol(contraction_flag_d, &contraction_flag, sizeof(int));
    cudaMemcpyToSymbol(contraction_prescription_d, &contraction_prescription, sizeof(int));
    
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
}
    
void CopyGlobalsToDevice2(void)
{
    cudaMemcpyToSymbol(potcor_d, &potcor, sizeof(double));
    cudaMemcpyToSymbol(psi_crit_d, &psi_crit, sizeof(double));
    cudaMemcpyToSymbol(psi_d_d, &psi_d, sizeof(double));
    cudaMemcpyToSymbol(psi_0_d, &psi_0, sizeof(double));
    cudaMemcpyToSymbol(log_rj2_d, &log_rj2, sizeof(double));
    cudaMemcpyToSymbol(psi_upper_initial_d, &psi_upper_initial, sizeof(double));
    cudaMemcpyToSymbol(psi_lower_initial_d, &psi_lower_initial, sizeof(double));
    cudaMemcpyToSymbol(Plcon_d, Plcon, 100*sizeof(double));
}

void DiskPotentialEstimateCUDA(void)
{
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> D_Dens_CUDA = thrust::device_malloc<double>(nr);
    
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    
    thrust::transform(Radius_CUDA, Radius_CUDA+nr, D_Dens_CUDA,
                      DiskDensityEstimateFunctor());
    
    thrust::copy(D_Dens_CUDA, D_Dens_CUDA+nr, D_Dens.begin());
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(D_Dens_CUDA);
}
 
void GasDiskCUDA(void)
{
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> GasDensity_Const_CUDA = thrust::device_malloc<double>(gas_disk*nr);
    thrust::device_ptr<double> Polytrope_Const_CUDA = thrust::device_malloc<double>(gas_disk*nr);
    
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    
    for (int j = 0; j < gas_disk; ++j)
    {
        int offset = j*(nr+1); 
        thrust::copy(Polytrope_Const[j].begin(), Polytrope_Const[j].end(), 
                     Polytrope_Const_CUDA+offset); 
        thrust::copy(GasDensity_Const[j].begin(), GasDensity_Const[j].end(),
                     GasDensity_Const_CUDA+offset); 
    }    
    
    for (int j = 0; j < gas_disk; ++j)
    {
        int offset = j*(nr+1); 
        thrust::transform(Radius_CUDA, Radius_CUDA+nr, Polytrope_Const_CUDA+offset,
                          GasDiskPolytropeFunctor(j, Radius_CUDA, A_Pot_CUDA, 
                                                  Polytrope_Const_CUDA, GasDensity_Const_CUDA));
        thrust::transform(Radius_CUDA, Radius_CUDA+nr, GasDensity_Const_CUDA+offset,
                          GasDiskDensityFunctor(j, Radius_CUDA, A_Pot_CUDA, 
                                                Polytrope_Const_CUDA, GasDensity_Const_CUDA));
    }    
    
    for (int j = 0; j < gas_disk; ++j)
    {
        int offset1 = j*(nr+1); 
        int offset2 = (j+1)*nr; 
        thrust::copy(Polytrope_Const_CUDA+offset1, Polytrope_Const_CUDA+offset2,
                     Polytrope_Const[j].begin()); 
        thrust::copy(GasDensity_Const_CUDA+offset1, GasDensity_Const_CUDA+offset2,
                     GasDensity_Const[j].begin()); 
    }    
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(Polytrope_Const_CUDA);
    thrust::device_free(GasDensity_Const_CUDA);
        
    for (int i = 0; i < nr; i+=1000)
    {
        cout << "Did " << i << " calculation " << setw(12) << Radius[i]
             << " " << Polytrope_Const[0][i] << "               "
             << Pot(Radius[i], 0) << " " << endl;
    }
    
    cout << "Gas Density:" << endl;
    
//     for (int i = 0; i < nr; i+=1000)
//     {
//         cout << "Did " << i << " calculation " << setw(12) << Radius[i]
//              << " " << GasDensity_Const_CUDA[0][i] << "               "
//              << Pot(Radius[i], 0) << " " << endl;
//     }
}

void GetADensCUDA(const int l)
{
    cudaMemcpyToSymbol(l_d, &l, sizeof(int));
    
    cudaMemcpyToSymbol(lmax_d, &lmax, sizeof(int));
    cudaMemcpyToSymbol(Q_total_d, &Q_total, sizeof(double));
    cudaMemcpyToSymbol(a00_d, &A_Pot[0][0], sizeof(double));
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> A_Dens_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Halo_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> GasDensity_Const_CUDA = thrust::device_malloc<double>(gas_disk*nr);
    thrust::device_ptr<double> Polytrope_Const_CUDA = thrust::device_malloc<double>(gas_disk*nr);

    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);
    thrust::copy(A_Dens[l/2].begin(), A_Dens[l/2].end(), A_Dens_CUDA);
    
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    
    thrust::copy(Dens_Psi_Halo.begin(), Dens_Psi_Halo.end(), Dens_Psi_Halo_CUDA);
    thrust::copy(Dens_Psi_Bulge.begin(), Dens_Psi_Bulge.end(), Dens_Psi_Bulge_CUDA);

    for (int j = 0; j < gas_disk; ++j)
    {
        int offset = j*(nr+1); 
        thrust::copy(Polytrope_Const[j].begin(), Polytrope_Const[j].end(), 
                     Polytrope_Const_CUDA+offset); 
        thrust::copy(GasDensity_Const[j].begin(), GasDensity_Const[j].end(),
                     GasDensity_Const_CUDA+offset); 
    }    
    
    thrust::transform(Radius_CUDA, Radius_CUDA+nr, A_Dens_CUDA, 
                      TotalDensityFunctor(Radius_CUDA, A_Pot_CUDA,
                                          Dens_Psi_Halo_CUDA, Dens_Psi_Bulge_CUDA, 
                                          Polytrope_Const_CUDA, GasDensity_Const_CUDA));
    
    thrust::copy(A_Dens_CUDA, A_Dens_CUDA+nr, A_Dens[l/2].begin());
    
//     for (int j = 0; j < gas_disk; ++j)
//     {
//         int offset1 = j*(nr+1); 
//         int offset2 = (j+1)*nr; 
//         thrust::copy(Polytrope_Const_CUDA+offset1, Polytrope_Const_CUDA+offset2,
//                      Polytrope_Const[j].begin()); 
//         thrust::copy(GasDensity_Const_CUDA+offset1, GasDensity_Const_CUDA+offset2,
//                      GasDensity_Const[j].begin()); 
//     }    
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(A_Dens_CUDA);
    thrust::device_free(Dens_Psi_Halo_CUDA);
    thrust::device_free(Dens_Psi_Bulge_CUDA);
    thrust::device_free(Polytrope_Const_CUDA);
    thrust::device_free(GasDensity_Const_CUDA);
}

void CopyGlobalsToDeviceChiSquare(void)
{
    cudaMemcpyToSymbol(disk_params_d, &disk_params, sizeof(int));
    cudaMemcpyToSymbol(nondisk_params_d, &nondisk_params, sizeof(int));
    cudaMemcpyToSymbol(astro_params_d, &astro_params, sizeof(int));
    cudaMemcpyToSymbol(error_params_d, &error_params, sizeof(int));
    cudaMemcpyToSymbol(do_chisq_file_io_d, &do_chisq_file_io, sizeof(int));
    cudaMemcpyToSymbol(n_los_d, &n_los, sizeof(int));
    cudaMemcpyToSymbol(n_vel_d, &n_vel, sizeof(int));
    
    cudaMemcpyToSymbol(total_mass_d, &total_mass, sizeof(double));
    cudaMemcpyToSymbol(halo_mass_d, &halo_mass, sizeof(double));
    cudaMemcpyToSymbol(total_disk_mass_d, &total_disk_mass, sizeof(double));
    cudaMemcpyToSymbol(bulge_mass_d, &bulge_mass, sizeof(double));
    cudaMemcpyToSymbol(r_edge_d, &r_edge, sizeof(double));
    cudaMemcpyToSymbol(halo_edge_d, &halo_edge, sizeof(double));
    cudaMemcpyToSymbol(disk_edge_d, &disk_edge, sizeof(double));
    cudaMemcpyToSymbol(bulge_edge_d, &bulge_edge, sizeof(double));
    
    cudaMemcpyToSymbol(ML_Disk, &Astro.ML_Disk[0], sizeof(Astro.ML_Disk));
    cudaMemcpyToSymbol(ml_bulge_d, &Astro.ml_bulge, sizeof(double));
    cudaMemcpyToSymbol(inclination_d, &Astro.inclination, sizeof(double));
    cudaMemcpyToSymbol(distance_d, &Astro.distance, sizeof(double));
    cudaMemcpyToSymbol(v_sys_d, &Astro.v_sys, sizeof(double));
    cudaMemcpyToSymbol(r_sys_d, &Astro.r_sys, sizeof(double));
    
    cudaMemcpyToSymbol(sin_inclination_d, &sin_inclination, sizeof(double));
    cudaMemcpyToSymbol(cos_inclination_d, &cos_inclination, sizeof(double));
}

double GetSurfaceBrightnessCUDA(const double error_factor, vector<double>& Radii, 
                                vector<double>& Data, vector<double>& Error, 
                                vector<double>& Ellipticity)
{
    thrust::device_vector<double> Radii_d(Radii.size());
    thrust::device_vector<double> Data_d(Radii.size());
    thrust::device_vector<double> Error_d(Radii.size());
    thrust::device_vector<double> Ellipticity_d(Radii.size());
    thrust::device_vector<double> chi2(Radii.size());

    cout << "here cuda " << Radii.size() << " " << Radii_d[0] << endl;
    thrust::copy(Radii.begin(), Radii.end(), Radii_d.begin());
    thrust::copy(Data.begin(), Data.end(), Data_d.begin());
    thrust::copy(Error.begin(), Error.end(), Error_d.begin());
    thrust::copy(Ellipticity.begin(), Ellipticity.end(), Ellipticity_d.begin());
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> F_R_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> DF_Sersic_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA = thrust::device_malloc<double>(nr);
    
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);      
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(F_R[m/2].begin(), F_R[m/2].end(), F_R_CUDA+offset);
    }
    thrust::copy(DF_Sersic.begin(), DF_Sersic.end(), DF_Sersic_CUDA);      
    thrust::copy(Dens_Psi_Bulge.begin(), Dens_Psi_Bulge.end(), Dens_Psi_Bulge_CUDA);      
    
    //should be possible to use transform reduce on what's below. on TODO list
    thrust::transform(make_zip_iterator(make_tuple(Radii_d.begin(), Data_d.begin(), 
                                                   Error_d.begin(), Ellipticity_d.begin())),
                      make_zip_iterator(make_tuple(Radii_d.end(), Data_d.end(), 
                                                   Error_d.end(), Ellipticity_d.end())),
                      chi2.begin(),
                      SurfaceBrightnessFunctor(error_factor, Radius_CUDA, 
                                               A_Pot_CUDA, F_R_CUDA, DF_Sersic_CUDA, 
                                               Dens_Psi_Bulge_CUDA));
    
    double chi_square = thrust::reduce(chi2.begin(), chi2.end(), 0, thrust::plus<double>());
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(DF_Sersic_CUDA);
    thrust::device_free(Dens_Psi_Bulge_CUDA);
    
    return chi_square;
}

double GetStarVelocityCUDA(const double error_factor, vector<double>& Radii, 
                           vector<double>& Data, vector<double>& Error)
{
    ////////////////////////////////////////////////////////////////////////////
    //size_t *pValue0, *pValue1, *pValue2;
    //cout << "here cudat" << endl;
    //cudaDeviceGetLimit(pValue0, cudaLimitStackSize);
    //cout << "here cudat" << endl;
    //cudaDeviceGetLimit(pValue1, cudaLimitPrintfFifoSize);
    //cout << "here cudat" << endl;
    //cudaDeviceGetLimit(pValue2, cudaLimitMallocHeapSize);
    //cout << "size  " << *pValue0 << "  " << *pValue1 << "  " << *pValue2 << endl;
    ////////////////////////////////////////////////////////////////////////////
    
    //cudaDeviceSetLimit(cudaLimitMallocHeapSize, 500000000);

    thrust::device_vector<double> Radii_d(Radii.size());
    thrust::device_vector<double> Data_d(Radii.size());
    thrust::device_vector<double> Error_d(Radii.size());
    thrust::device_vector<double> chi2(Radii.size());

    cout << "here cuda" << endl;
    thrust::copy(Radii.begin(), Radii.end(), Radii_d.begin());
    thrust::copy(Data.begin(), Data.end(), Data_d.begin());
    thrust::copy(Error.begin(), Error.end(), Error_d.begin());
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> F_R_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> Omega_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Omega2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_K_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_K2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> DF_Sersic_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Am_Tab_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab2_Zero_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Rad_Spline_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> FD_Rat_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> D_Rat2_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> FSZ_Rat_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> SZ_Rat2_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);  
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(F_R[m/2].begin(), F_R[m/2].end(), F_R_CUDA+offset);
    }
    thrust::copy(Omega.begin(), Omega.end(), Omega_CUDA);      
    thrust::copy(Omega2.begin(), Omega2.end(), Omega2_CUDA);      
    thrust::copy(A_K.begin(), A_K.end(), A_K_CUDA);      
    thrust::copy(A_K2.begin(), A_K2.end(), A_K2_CUDA);      
    thrust::copy(DF_Sersic.begin(), DF_Sersic.end(), DF_Sersic_CUDA);      
    thrust::copy(Dens_Psi_Bulge.begin(), Dens_Psi_Bulge.end(), Dens_Psi_Bulge_CUDA);      
    thrust::copy(Am_Tab.begin(), Am_Tab.end(), Am_Tab_CUDA);      
    thrust::copy(R_Tab.begin(), R_Tab.end(), R_Tab_CUDA);      
    thrust::copy(R_Tab2.begin(), R_Tab2.end(), R_Tab2_CUDA);      
    thrust::copy(R_Tab2_Zero.begin(), R_Tab2_Zero.end(), R_Tab2_Zero_CUDA);      
    cout << "here cuda1" << endl;
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1); 
        thrust::copy(Rad_Spline[j].begin(), Rad_Spline[j].end(), Rad_Spline_CUDA+offset); 
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(FD_Rat[j].begin(), FD_Rat[j].end(), FD_Rat_CUDA+offset);      
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(D_Rat2[j].begin(), D_Rat2[j].end(), D_Rat2_CUDA+offset);      
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(FSZ_Rat[j].begin(), FSZ_Rat[j].end(), FSZ_Rat_CUDA+offset);      
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(SZ_Rat2[j].begin(), SZ_Rat2[j].end(), SZ_Rat2_CUDA+offset);      
    }    
    
    //should be possible to use transform reduce on what's below. on TODO list
    thrust::transform(make_zip_iterator(make_tuple(Radii_d.begin(), Data_d.begin(), Error_d.begin())),
                      make_zip_iterator(make_tuple(Radii_d.end(), Data_d.end(), Error_d.end())),
                      chi2.begin(),
                      StarVelocityFunctor(error_factor, Radius_CUDA, Rad_Spline_CUDA,
                                          A_Pot_CUDA, F_R_CUDA, Omega_CUDA, Omega2_CUDA, A_K_CUDA, 
                                          A_K2_CUDA, DF_Sersic_CUDA, Dens_Psi_Bulge_CUDA, 
                                          Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_CUDA, R_Tab2_Zero_CUDA,  
                                          FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA));
    
    double chi_square = thrust::reduce(chi2.begin(), chi2.end(), 0, thrust::plus<double>());
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(F_R_CUDA);
    thrust::device_free(Omega_CUDA);
    thrust::device_free(Omega2_CUDA);
    thrust::device_free(A_K_CUDA);
    thrust::device_free(A_K2_CUDA);
    thrust::device_free(DF_Sersic_CUDA);
    thrust::device_free(Dens_Psi_Bulge_CUDA);
    thrust::device_free(Am_Tab_CUDA);
    thrust::device_free(R_Tab_CUDA);
    thrust::device_free(R_Tab2_CUDA);
    thrust::device_free(R_Tab2_Zero_CUDA);
    thrust::device_free(Rad_Spline_CUDA);
    thrust::device_free(FD_Rat_CUDA);
    thrust::device_free(D_Rat2_CUDA);
    thrust::device_free(FSZ_Rat_CUDA);
    thrust::device_free(SZ_Rat2_CUDA);
    
    return chi_square;
}

double GetStarDispersionCUDA(const double error_factor, vector<double>& Radii, 
                             vector<double>& Data, vector<double>& Error)
{
    thrust::device_vector<double> Radii_d(Radii.size());
    thrust::device_vector<double> Data_d(Radii.size());
    thrust::device_vector<double> Error_d(Radii.size());
    thrust::device_vector<double> chi2(Radii.size());

    thrust::copy(Radii.begin(), Radii.end(), Radii_d.begin());
    thrust::copy(Data.begin(), Data.end(), Data_d.begin());
    thrust::copy(Error.begin(), Error.end(), Error_d.begin());
    
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> F_R_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> Omega_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Omega2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_K_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_K2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> DF_Sersic_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Dens_Psi_Bulge_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Am_Tab_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab2_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> R_Tab2_Zero_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> Rad_Spline_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> FD_Rat_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> D_Rat2_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> FSZ_Rat_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    thrust::device_ptr<double> SZ_Rat2_CUDA = thrust::device_malloc<double>(disk*(nr_spline+1));
    
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);      
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(F_R[m/2].begin(), F_R[m/2].end(), F_R_CUDA+offset);
    }
    thrust::copy(Omega.begin(), Omega.end(), Omega_CUDA);      
    thrust::copy(Omega2.begin(), Omega2.end(), Omega2_CUDA);      
    thrust::copy(A_K.begin(), A_K.end(), A_K_CUDA);      
    thrust::copy(A_K2.begin(), A_K2.end(), A_K2_CUDA);      
    thrust::copy(DF_Sersic.begin(), DF_Sersic.end(), DF_Sersic_CUDA);      
    thrust::copy(Dens_Psi_Bulge.begin(), Dens_Psi_Bulge.end(), Dens_Psi_Bulge_CUDA);      
    thrust::copy(Am_Tab.begin(), Am_Tab.end(), Am_Tab_CUDA);      
    thrust::copy(R_Tab.begin(), R_Tab.end(), R_Tab_CUDA);      
    thrust::copy(R_Tab2.begin(), R_Tab2.end(), R_Tab2_CUDA);      
    thrust::copy(R_Tab2_Zero.begin(), R_Tab2_Zero.end(), R_Tab2_Zero_CUDA);      
    //thrust::copy(Rad_Spline[j].begin(), Rad_Spline[j].end(), Rad_Spline_CUDA);      
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(Rad_Spline[j].begin(), Rad_Spline[j].end(), Rad_Spline_CUDA+offset);  
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(FD_Rat[j].begin(), FD_Rat[j].end(), FD_Rat_CUDA+offset);  
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(D_Rat2[j].begin(), D_Rat2[j].end(), D_Rat2_CUDA+offset);      
    }    
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(FSZ_Rat[j].begin(), FSZ_Rat[j].end(), FSZ_Rat_CUDA+offset);
    }     
    for (int j = 0; j < disk; ++j)
    {
        int offset = j*(nr_spline+1);
        thrust::copy(SZ_Rat2[j].begin(), SZ_Rat2[j].end(), SZ_Rat2_CUDA+offset);   
    }   
    
    //should be possible to use transform reduce on what's below. on TODO list
    thrust::transform(make_zip_iterator(make_tuple(Radii_d.begin(), Data_d.begin(), Error_d.begin())),
                      make_zip_iterator(make_tuple(Radii_d.end(), Data_d.end(), Error_d.end())),
                      chi2.begin(),
                      StarDispersionFunctor(error_factor, Radius_CUDA, Rad_Spline_CUDA, 
                                            A_Pot_CUDA, F_R_CUDA, Omega_CUDA, 
                                            Omega2_CUDA, A_K_CUDA, A_K2_CUDA, DF_Sersic_CUDA,
                                            Dens_Psi_Bulge_CUDA, Am_Tab_CUDA,
                                            R_Tab_CUDA, R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                                            FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA));
    
    double chi_square = thrust::reduce(chi2.begin(), chi2.end(), 0, thrust::plus<double>());
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(F_R_CUDA);
    thrust::device_free(Omega_CUDA);
    thrust::device_free(Omega2_CUDA);
    thrust::device_free(A_K_CUDA);
    thrust::device_free(A_K2_CUDA);
    thrust::device_free(DF_Sersic_CUDA);
    thrust::device_free(Dens_Psi_Bulge_CUDA);
    thrust::device_free(Am_Tab_CUDA);
    thrust::device_free(R_Tab_CUDA);
    thrust::device_free(R_Tab2_CUDA);
    thrust::device_free(R_Tab2_Zero_CUDA);
    thrust::device_free(Rad_Spline_CUDA);
    thrust::device_free(FD_Rat_CUDA);
    thrust::device_free(D_Rat2_CUDA);
    thrust::device_free(FSZ_Rat_CUDA);
    thrust::device_free(SZ_Rat2_CUDA);
    
    return chi_square;
}

double GetCircVelocityCUDA(const double error_factor, vector<double>& Radii, 
                           vector<double>& Data, vector<double>& Error)
{
    thrust::device_vector<double> Radii_d(Radii.size());
    thrust::device_vector<double> Data_d(Radii.size());
    thrust::device_vector<double> Error_d(Radii.size());
    thrust::device_vector<double> chi2(Radii.size());

    thrust::copy(Radii.begin(), Radii.end(), Radii_d.begin());
    thrust::copy(Data.begin(), Data.end(), Data_d.begin());
    thrust::copy(Error.begin(), Error.end(), Error_d.begin());
    
    //cout << "here cuda" << endl;
    thrust::device_ptr<double> Radius_CUDA = thrust::device_malloc<double>(nr);
    thrust::device_ptr<double> A_Pot_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
    thrust::device_ptr<double> F_R_CUDA = thrust::device_malloc<double>(nr*l_max/2+nr);
      
    thrust::copy(Radius.begin(), Radius.end(), Radius_CUDA);      
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(A_Pot[m/2].begin(), A_Pot[m/2].end(), A_Pot_CUDA+offset);
    }
    for (int m = 0; m < l_max+1; m+=2)
    {
        int offset = m*nr/2;
        thrust::copy(F_R[m/2].begin(), F_R[m/2].end(), F_R_CUDA+offset);
    }
    
    //should be possible to use transform reduce on what's below. on TODO list
    thrust::transform(make_zip_iterator(make_tuple(Radii_d.begin(), Data_d.begin(), Error_d.begin())),
                      make_zip_iterator(make_tuple(Radii_d.end(), Data_d.end(), Error_d.end())),
                      chi2.begin(),
                      CircVelocityFunctor(error_factor, Radius_CUDA, A_Pot_CUDA, F_R_CUDA));
    
    //double chi_square = thrust::reduce(chi2.begin(), chi2.end(), 0, plus<double>);
    
    thrust::device_free(Radius_CUDA);
    thrust::device_free(A_Pot_CUDA);
    thrust::device_free(F_R_CUDA);
    
    return thrust::reduce(chi2.begin(), chi2.end(), 0, thrust::plus<double>());
    //return chi2[0];
}

// struct RandomNumber
// {
//     __host__ __device__
//     float operator()(int index)
//     {
//         default_random_engine rng;
// 
//         rng.discard(2*index);
// 
//         return (float)rng()/default_random_engine::max;
//     }
// };

