#include "galaxy_cuda.h"

__host__ __device__   
double HaloDensityCUDA(double &r, 
                       thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                       thrust::device_ptr<double> Halo_AC_Dens_CUDA)
{
    double nfwdens = HaloProfileDensCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA);
    
    double trunc_fac = GetHaloTruncCUDA(r);
    //if(r<0.01)printf("denshalo %.12f %.12f %.12f\n", r, nfwdens, trunc_fac);
    
    return nfwdens*trunc_fac;
}

__host__ __device__   
double HaloDensityPrimeCUDA(double &r, 
                            thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                            thrust::device_ptr<double> Halo_AC_Dens_CUDA, 
                            thrust::device_ptr<double> Halo_AC_Dens_D_CUDA)
{
    double nfwdens = HaloProfileDensCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA);
    double dnfwdens = HaloProfileDensPrimeCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_D_CUDA);
    
    double trunc_fac = GetHaloTruncCUDA(r);
    double dtrunc_fac = GetHaloTruncPrimeCUDA(r);
    
    return nfwdens*dtrunc_fac + dnfwdens*trunc_fac;
}

__host__ __device__   
double HaloDensity2PrimeCUDA(double &r, 
                             thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_D_CUDA, 
                             thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA)
{
    double nfwdens = HaloProfileDensCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_CUDA);
    double dnfwdens = HaloProfileDensPrimeCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_D_CUDA);
    double ddnfwdens = HaloProfileDens2PrimeCUDA(r, Halo_AC_Radius_CUDA, Halo_AC_Dens_DD_CUDA);
    
    double trunc_fac = GetHaloTruncCUDA(r);
    double dtrunc_fac = GetHaloTruncPrimeCUDA(r);
    double ddtrunc_fac = GetHaloTrunc2PrimeCUDA(r);
    
    //cout << "h2prime  " << nfwdens << " " << dnfwdens << " " << ddnfwdens << " "
    //     << trunc_fac << " " << dtrunc_fac << " " << ddtrunc_fac << endl;
    //cout <<  ddnfwdens*trunc_fac + 2*dnfwdens*dtrunc_fac + nfwdens*ddtrunc_fac << endl;
    
    return ddnfwdens*trunc_fac + 2*dnfwdens*dtrunc_fac + nfwdens*ddtrunc_fac;
}
    
__host__ __device__   
double RawHaloProfileCUDA(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/a_halo_d;
    
    double density = halo_const_d / pow(s,cusp_d) / pow(1+s,3-cusp_d);
    
    //if(r<0.01)printf("denshalo %.12f %.12f %.12f\n", r, halo_const_d, density);
    return density;
    
//     double s = r/a_halo_d;
//     double arg = -2.0*(pow(s, cusp_d)-1.0)/cusp_d;
//     
//     double density = v_halo_d*exp(arg);
//     
//     //cout << "          " << r << " " << density << endl;
//     
//     return density;
}

__host__ __device__   
double RawHaloProfilePrimeCUDA(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/a_halo_d;
    
    double ddensity = -RawHaloProfileCUDA(r)*(3*s + cusp_d)/a_halo_d/s/(1 + s);
    
    return ddensity;
    
//     double s = r/a_halo_d;
//     double arg = -2.0*(pow(s, cusp_d)-1.0)/cusp_d;
//     
//     double ddensity = -2*v_halo_d*exp(arg)*pow(r, cusp_d-1)/pow(a_halo_d, cusp_d);
//     
//     return ddensity;
}

__host__ __device__   
double RawHaloProfile2PrimeCUDA(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/a_halo_d;
    
    double deno = s*s*a_halo_d*a_halo_d*(1+s)*(1+s);
    double dddensity = RawHaloProfileCUDA(r)*(cusp_d*(cusp_d+1) + 8*cusp_d*s + 12*s*s) / deno;
    
    return dddensity;
    
//     double s = r/a_halo_d;
//     double arg = -2.0*(pow(s, cusp_d)-1.0)/cusp_d;
//     double exparg = exp(arg);
//     double recip = 1.0/pow(a_halo_d, cusp_d);
//     double temp1 = (cusp_d-1.0)*pow(r, cusp_d-2.0);
//     double temp2 = 2.0*recip*pow(r, 2.0*cusp_d-2.0);
//     
//     double dddensity = -2*v_halo_d*exparg*recip*(temp1-temp2);
//     
//     return dddensity;
}

__host__ __device__   
double GetHaloTruncCUDA(double &r)
{
    double truncfac;
    double erfarg = (r - c_halo_d) * oneoversqrt2 / drtrunc_halo_d;
    
    if (erfarg > 6) return truncfac = 0;
    else if (erfarg < -6) return truncfac = 1;
    else return truncfac = 0.5 * erfc(erfarg);
}

__host__ __device__   
double GetHaloTruncPrimeCUDA(double &r)
{
    double truncfac;
    double t = (r - c_halo_d) * oneoversqrt2 / drtrunc_halo_d;
    t *= t;
    
    if (t > 36) return truncfac = 0;
    else return truncfac = -exp(-t)/drtrunc_halo_d*oneoversqrt2pi;
}

__host__ __device__   
double GetHaloTrunc2PrimeCUDA(double &r)
{
    double truncfac;
    double t = (r - c_halo_d) * oneoversqrt2 / drtrunc_halo_d;
    double tt = t*t;
    
    if (tt > 36) return truncfac = 0;
    else return truncfac = t*exp(-tt)*oneoversqrtpi/drtrunc_halo_d/drtrunc_halo_d;
}

__host__ __device__   
double HaloForceCUDA(double &r, thrust::device_ptr<double> Radius_CUDA, 
                     thrust::device_ptr<double> H_FR_CUDA)
{
    double log_r = log10(r);
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
    
    return t*H_FR_CUDA[ihi] + tm1*H_FR_CUDA[ihi-1];
}

__host__ __device__
double HaloProfileDensCUDA(double &r, thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                           thrust::device_ptr<double> Halo_AC_Dens_CUDA)
{
    if (!contraction_flag_d)
    {
        return RawHaloProfileCUDA(r);
    }
    else if (contraction_flag_d)
    {
        int rad_index = 0;
        
        if (r <= Halo_AC_Radius_CUDA[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius_CUDA[nr_ac_d-1])
        {
            return RawHaloProfileCUDA(r);//*(1-baryon_frac);
        }
        else
        {
            while (Halo_AC_Radius_CUDA[rad_index] < r)
            {
                ++rad_index;
            }
        }
        
        if (rad_index < 1)
        {
            printf("rad_index is < 1 in HaloProfileDens! Everything after this");
            printf(" point is probably garbage.\n");
        }
        else if (rad_index > nr_ac_d-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac_d-1;
        }

        double r1 = Halo_AC_Radius_CUDA[rad_index-1];
        double r2 = Halo_AC_Radius_CUDA[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << "Ac " << rad_index << " " << r << " " << t << " " << Halo_AC_Dens[rad_index]
        //     << "     " << tm1 << " " << Halo_AC_Dens[rad_index-1] << endl;
        
        return t*Halo_AC_Dens_CUDA[rad_index] + tm1*Halo_AC_Dens_CUDA[rad_index-1];  
         
//     int ihi = ceil(r/dr_d);
//     
//     if(r < dr_d)
//     {
//         ihi = 1;
//     }
//     else if (ihi < 1)
//     {
//     }
//     else if (ihi > nr_ac_d-1)
//     {
//         //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
//         //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
//         ihi = nr_d-1;
//     }
//     
//     double r1 = Halo_AC_Radius_CUDA[ihi-1];
//     double r2 = Halo_AC_Radius_CUDA[ihi];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     return t*Halo_AC_Dens_CUDA[ihi] + tm1*Halo_AC_Dens_CUDA[ihi-1];
    }
}  
        
__host__ __device__
double HaloProfileDensPrimeCUDA(double &r, thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                                thrust::device_ptr<double> Halo_AC_Dens_D_CUDA)
{
    if (!contraction_flag_d)
    {
        return RawHaloProfilePrimeCUDA(r);
    }
    else if (contraction_flag_d)
    {
        int rad_index = 0;
        
        if (r <= Halo_AC_Radius_CUDA[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius_CUDA[nr_ac_d-1])
        {
            return RawHaloProfilePrimeCUDA(r);//*(1-baryon_frac);
        }
        else
        {
            while (Halo_AC_Radius_CUDA[rad_index] < r)
            {
                ++rad_index;
            }
        }
        
        if (rad_index < 1)
        {
            printf("rad_index is < 1 in HaloProfileDensPrime! Everything after this");
            printf(" point is probably garbage.\n");
        }
        else if (rad_index > nr_ac_d-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac_d-1;
        }

        double r1 = Halo_AC_Radius_CUDA[rad_index-1];
        double r2 = Halo_AC_Radius_CUDA[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << "Ac " << rad_index << " " << r << " " << t << " " << Halo_AC_Dens[rad_index]
        //     << "     " << tm1 << " " << Halo_AC_Dens[rad_index-1] << endl;
        
        return t*Halo_AC_Dens_D_CUDA[rad_index] + tm1*Halo_AC_Dens_D_CUDA[rad_index-1];   

//     int ihi = ceil(r/dr_d);
//     
//     if(r < dr_d)
//     {
//         ihi = 1;
//     }
//     else if (ihi < 1)
//     {
//     }
//     else if (ihi > nr_ac_d-1)
//     {
//         //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
//         //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
//         ihi = nr_d-1;
//     }
//     
//     double r1 = Halo_AC_Radius_CUDA[ihi-1];
//     double r2 = Halo_AC_Radius_CUDA[ihi];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     return t*Halo_AC_Dens_D_CUDA[ihi] + tm1*Halo_AC_Dens_D_CUDA[ihi-1];
    }
}  
    
__host__ __device__
double HaloProfileDens2PrimeCUDA(double &r, thrust::device_ptr<double> Halo_AC_Radius_CUDA, 
                                 thrust::device_ptr<double> Halo_AC_Dens_DD_CUDA)
{
    if (!contraction_flag_d)
    {
        return RawHaloProfile2PrimeCUDA(r);
    }
    else if (contraction_flag_d)
    {
        int rad_index = 0;
        
        if (r <= Halo_AC_Radius_CUDA[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius_CUDA[nr_ac_d-1])
        {
            return RawHaloProfile2PrimeCUDA(r);//*(1-baryon_frac);
        }
        else
        {
            while (Halo_AC_Radius_CUDA[rad_index] < r)
            {
                ++rad_index;
            }
        }
        
        if (rad_index < 1)
        {
            printf("rad_index is < 1 in HaloProfileDens2Prime! Everything after this");
            printf(" point is probably garbage.\n");
        }
        else if (rad_index > nr_ac_d-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac_d-1;
        }

        double r1 = Halo_AC_Radius_CUDA[rad_index-1];
        double r2 = Halo_AC_Radius_CUDA[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << "Ac " << rad_index << " " << r << " " << t << " " << Halo_AC_Dens[rad_index]
        //     << "     " << tm1 << " " << Halo_AC_Dens[rad_index-1] << endl;
        
        return t*Halo_AC_Dens_DD_CUDA[rad_index] + tm1*Halo_AC_Dens_DD_CUDA[rad_index-1]; 
  
//     int ihi = ceil(r/dr_d);
//     
//     if(r < dr_d)
//     {
//         ihi = 1;
//     }
//     else if (ihi < 1)
//     {
//     }
//     else if (ihi > nr_ac_d-1)
//     {
//         //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
//         //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
//         ihi = nr_d-1;
//     }
//     
//     double r1 = Halo_AC_Radius_CUDA[ihi-1];
//     double r2 = Halo_AC_Radius_CUDA[ihi];
//     double t = (r-r1)/(r2-r1);
//     double tm1 = 1-t;
//     
//     return t*Halo_AC_Dens_DD_CUDA[ihi] + tm1*Halo_AC_Dens_DD_CUDA[ihi-1];
    }
}  
