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
    //printf("Nowpot %f\n", r);
    double r1 = *(Radius_CUDA+ihi-1);
    double r2 = *(Radius_CUDA+ihi);
    //printf("Nowpot1 %f\n", r);
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    double cos_theta = z/r;
    
    for (int l = l_max_d; l > -1; l-=2)
    {
        pot += LegendreCUDA(l, cos_theta)*Plcon_d[l]*
               (t*(*(A_Pot_CUDA+l/2*nr_d+ihi)) + tm1*(*(A_Pot_CUDA+l/2*nr_d+ihi-1)));
    }
    
    if (disk_flag_d)
    {
        pot += AppDiskPotCUDA(s,z);
    }
        
    if (gasdisk_flag_d)
    {
        pot += AppGasDiskPotCUDA(s,z);
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

__host__ __device__     
double GetTruncGasCUDA(const double &rad, const int &i)
{
    //double truncfac;
    double erfarg = (rad - (*(Out_GasDisk+i))) * oneoversqrt2 / (*(Dr_Trunc_Gas+i));
    
    if (erfarg > 4) return 0;
    else if (erfarg < -4) return 1;
    else return 0.5 * erfc(erfarg);
}

__host__ __device__     
void GetTruncGasPrimeCUDA(const double &r, const int &i, 
                          double &truncfac, double &truncfacprime)
{
    double erfarg = (r - (*(Out_GasDisk+i)))*oneoversqrt2/(*(Dr_Trunc_Gas+i));
    
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
        truncfacprime = -exp(-erfarg*erfarg)/(*(Dr_Trunc_Gas+i))*oneoversqrt2pi;
    }
}

__host__ __device__     
double BulgeDensPsiCUDA(double energy, thrust::device_ptr<double> Dens_Psi_Bulge_CUDA)
{
    if (energy < psi_crit_d) return 0;
    
    else if (energy >= psi_0_d) return Dens_Psi_Bulge_CUDA[0];
    
    double rj1 = (psi_0_d - energy)/(psi_0_d - psi_d_d);
    //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
    double rj = (n_psi_d-1.0)*log(rj1)*log_rj2_d;
    
	int j = int(rj);
    
    if (j < 0)
        j = 0;
    else if (j >= n_psi_d-1)
        j = n_psi_d - 2;
    
    double frac = rj - j;
    
    //cout << "bdpsi " << rj1 << " " << rj2 << " " << rj << " " << j << 
    //        " " << setprecision(12) << psi_0 << " " << energy << endl;
    
    return Dens_Psi_Bulge_CUDA[j] + frac*(Dens_Psi_Bulge_CUDA[j+1] - Dens_Psi_Bulge_CUDA[j]);
}
    
__host__ __device__     
void ForceCUDA(double s, double z, double &force_s, double &force_z,
               thrust::device_ptr<double> Radius_CUDA,                            
               thrust::device_ptr<double> A_Pot_CUDA, 
               thrust::device_ptr<double> F_R_CUDA)
{ 
    double pc[20], p[20], dp[20];//, pot;
    
    double r=sqrt(s*s+z*z);
    
    int ihi = ceil(r/dr_d);
    int ihim = ihi - 1;
    
    if (r<dr_d)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        //There really should be an exit here. Must find out how to do that in CUDA.
    }
    else if (ihi > nr_d-1)
    {
        ihi = nr_d-1;
    }
    
    double r1 = *(Radius_CUDA+ihi-1);//Radius[ihi-1];
    double r2 = *(Radius_CUDA+ihi);//Radius[ihi];
    double r_edge = *(Radius_CUDA+nr_d-1);
    double t = (r-r1)/(r2-r1);
    double tm = 1-t;
    
    if (r==0)
    {
        force_s = force_z = 0;
    }
    else
    {
        double cos_theta = z/r;
        double cos2theta = cos_theta*cos_theta;
        double sin_theta = s/r;
        double sin2theta = 1 - cos2theta;
        
        for (int l = 0; l < lmax_d+1; l+=2)
        {
            pc[l/2] = sqrt((2.0*l+1)*oneover4pi);
            p[l/2] = LegendreCUDA(l, cos_theta);
            
            if (fabs(cos_theta) == 1)
            {
                dp[l/2] = 0;
            }
            else
            {
                dp[l/2] = l*(LegendreCUDA(l, cos_theta)-cos_theta*p[l/2])/sin2theta;
            }   
        }
        
        for (int i = 0; i < lmax_d/2 + 1; ++i)
        {
            p[i] = p[i]*pc[i];
            dp[i] = dp[i]*pc[i];
        }
        
        double frr = 0, fth = 0, pot = 0;
        
        if (r <= r_edge)
        {
            for (int i = 0; i < lmax_d/2 + 1; ++i)
            {
                frr += p[i]*(t*F_R_CUDA[i*nr_d+ihi] + tm*F_R_CUDA[i*nr_d+ihim]);
            }
            
            for (int i = 2; i < lmax_d/2 + 2; ++i)
            {
                fth -= sin_theta*dp[i-1]*(t*A_Pot_CUDA[(i-1)*nr_d+ihi] + tm*A_Pot_CUDA[(i-1)*nr_d+ihim]);
            }
            
            for (int i = 0; i < lmax_d/2 + 1; ++i)
            {
                pot += p[i]*(t*A_Pot_CUDA[i*nr_d+ihi] + tm*A_Pot_CUDA[i*nr_d+ihim]);
            }
            
        }
        else
        {
            for (int i = 0; i < lmax_d/2 + 1; ++i)
            {
                int l = 2*i;
                frr -= (l+1)*p[i]*A_Pot_CUDA[i*nr_d+nr_d] / r_edge * pow(r_edge/r,l+2);
            }
            
            for (int i = 2; i < lmax_d/2 + 2; ++i)
            {
                int l = 2*(i-1);
                fth -= sin_theta*dp[i-1]*A_Pot_CUDA[(i-1)*nr_d+nr_d] * pow(r_edge/r,l+1);
            }
            
            for (int i = 0; i < lmax_d/2 + 1; ++i)
            {
                int l = 2*i;
                pot += p[i]*A_Pot_CUDA[i*nr_d+nr_d] * pow(r_edge/r,l+1);
            }
        }
        
        if (disk_flag_d == 1)
        {
            pot += AppDiskPotCUDA(s,z);
        }
        
        force_s = -sin_theta*frr - cos_theta*fth/r;
        force_z = -cos_theta*frr + sin_theta*fth/r;
        
        if (disk_flag_d == 1)
        {
            double fsad, fzad;
            AppDiskForceCUDA(s,z,fsad,fzad);
            force_s += fsad;
            force_z += fzad;
        }
    }
}    
        
__host__ __device__     
double BulgeDFCUDA(double &energy, thrust::device_ptr<double> DF_Sersic_CUDA)
{
    if (energy < psi_crit_d) return 0;
    
    if (energy >= psi_0_d) return exp(*(DF_Sersic_CUDA));
    
    double rj1 = (psi_0_d - energy)/(psi_0_d - psi_d_d);
    //double rj2 = (psi_0 - psi_crit)/(psi_0 - psi_d);
    double rj = 1 + (n_psi_d-1.0)*log(rj1)*log_rj2_d;// /log(rj2);
    
    int j = int(rj);
    
    if (j < 1)
        j = 1;
    else if (j >= n_psi_d)
        j = n_psi_d - 1;
    
    double frac = rj - j;
    
    return exp((*(DF_Sersic_CUDA+j-1)) + frac*((*(DF_Sersic_CUDA+j)) - (*(DF_Sersic_CUDA+j-1))));
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// DF functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

__host__ __device__     
void GetOmegaKappaCUDA(double &r, double &freq_omega, double &freq_kappa,
                       thrust::device_ptr<double> Radius_CUDA,
                       thrust::device_ptr<double> Omega_CUDA,
                       thrust::device_ptr<double> Omega2_CUDA,
                       thrust::device_ptr<double> A_K_CUDA,
                       thrust::device_ptr<double> A_K2_CUDA)
{
    // KK: a anrd b are coefficients of potential -a/r -b/r^3 +c 
    // which give omeka and kappa as observed
    // at the outer radius of the model---these are used to extrapolate if needed.
    // in this potential om^2=a/r^3+3b/r^5, kap^2=a/r^3-3b/r^5 
    double rad = Radius_CUDA[nr_d-1];
    double a = rad*rad*rad/2*(Omega_CUDA[nr_d-1]*Omega_CUDA[nr_d-1]+
                              A_K_CUDA[nr_d-1]*A_K_CUDA[nr_d-1]);
    double b = pow(rad, 5)/6*(Omega_CUDA[nr_d-1]*Omega_CUDA[nr_d-1]-
                              A_K_CUDA[nr_d-1]*A_K_CUDA[nr_d-1]);
    //double c = a/rad/rad/rad+b/pow(rad, 5)+Pot_Major_Tot[nr_d-1];
    
    if (r > Radius_CUDA[nr_d-1])
    {
        double r3 = a/rad/rad/rad;
        double r5 = 3*b/pow(rad, 5);
        freq_omega = sqrt(r3+r5);
        freq_kappa = sqrt(r3-r5);
    }
    else
    {
        SplintDCUDA(Radius_CUDA, Omega_CUDA, Omega2_CUDA, nr_d, r, freq_omega);
        SplintDCUDA(Radius_CUDA, A_K_CUDA, A_K2_CUDA, nr_d, r, freq_kappa);
    }
}

__host__ __device__     
void SplintDCUDA(const thrust::device_ptr<double> xa, const thrust::device_ptr<double> ya, 
                 const thrust::device_ptr<double> y2a, int n, double x, double &y)
{
    int k, khi = n, klo = 1;
	double h;
    
//     if (khi < n)
//     {
//         if (xa.at(khi)>x && xa.at(klo)<x && khi-klo==1)
//         {
//             h = xa.at(khi)-xa.at(klo);
//             if (h < 0)
//             {
//                 cerr << "bad xa input in splint. h < 0" << endl;
//                 return;
//             }
//         }
//     }
            
    klo = 0;
    khi = n-1;
            
    while (khi-klo > 1)
    {
        k = (khi+klo)/2;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }
        
    h = xa[khi]-xa[klo];
	
    if (h <= 0)
    {
        printf("Bad xa input in splint");//. h <= 0. %f %f %f\n", xa[khi], xa[klo], khi);
        return;
    }
        
    double a = (xa[khi]-x)/h;
    double b = (x-xa[klo])/h;
      
    y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*h*h/6;
        
    return;
}

__host__ __device__     
double DiskDF5ezCUDA(double &vr, double &vt, double &vz, double &r, double &z, int &i,
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
    double psi_r0 = PotCUDA(r, 0, Radius_CUDA, A_Pot_CUDA), psi_rz;
    
    if (z == 0)
    {
        psi_rz = psi_r0;
    }
    else
    {
        psi_rz = PotCUDA(r, z, Radius_CUDA, A_Pot_CUDA);
    }
    
    double ep = 0.5*(vr*vr+vt*vt) - psi_r0;
    double am = r*vt;
    double ez = 0.5*vz*vz - psi_rz + psi_r0;
    
    double d = DiskDF3ezCUDA(ep, am, ez, i, Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA, Omega_CUDA,
                             Omega2_CUDA, A_K_CUDA, A_K2_CUDA, Am_Tab_CUDA, R_Tab_CUDA, 
                             R_Tab2_CUDA, R_Tab2_Zero_CUDA, 
                             FD_Rat_CUDA, D_Rat2_CUDA, FSZ_Rat_CUDA, SZ_Rat2_CUDA);
    
    //cout << "df   " << d << " " << ep << " " << am << " " << ez << endl;
    
    return d;
}

__host__ __device__     
double DiskDF3ezCUDA(double &ep, double &am, double &ez, int &i,
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
    double psi00 = PotCUDA(0, 0, Radius_CUDA, A_Pot_CUDA);
    double rc = RCircCUDA(am, Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_Zero_CUDA);
    double freq_omega, freq_kappa;
    
    GetOmegaKappaCUDA(rc, freq_omega, freq_kappa, Radius_CUDA, Omega_CUDA, 
                      Omega2_CUDA, A_K_CUDA, A_K2_CUDA);
    
    double v_circ = rc*freq_omega;
    
    double ec = -PotCUDA(rc, 0, Radius_CUDA, A_Pot_CUDA) + 0.5*v_circ*v_circ;
    
    if (am < 0)
    {
        ec = -2*psi00-ec;
    }
    
    double sr2 = SigR2CUDA(rc,i);
    double sz2 = SigZ2CUDA(rc,i, Radius_CUDA, Rad_Spline_CUDA, A_Pot_CUDA,
                           FSZ_Rat_CUDA,SZ_Rat2_CUDA);
    
    if (sz2 > 0)
    {
        double f_vert = FnaMidDenCUDA(rc,i,Radius_CUDA,Rad_Spline_CUDA,A_Pot_CUDA,FD_Rat_CUDA,D_Rat2_CUDA)*
                        exp(-ez/sz2)*oneoversqrt2pi/sqrt(sz2);
        f_vert = f_vert*freq_omega/freq_kappa/sr2*exp((ec-ep)/sr2)*oneoverpi;
        
        //if(f_vert<0)cout<<"fvert<0  "<<freq_omega/freq_kappa/sr2<<" "<<FnaMidDenCUDA(rc,i)<<endl;
        
        return f_vert;
    }
    else
    {
        return 0;
    }
}
    
__host__ __device__     
double FnaMidDenCUDA(double &r, int &j, thrust::device_ptr<double> Radius_CUDA, 
                     thrust::device_ptr<double> Rad_Spline_CUDA,
                     thrust::device_ptr<double> A_Pot_CUDA, 
                     thrust::device_ptr<double> FD_Rat_CUDA,
                     thrust::device_ptr<double> D_Rat2_CUDA)
{
    double f_cor;
    
    SplintDCUDA(Rad_Spline_CUDA+j*nr_spline_d, FD_Rat_CUDA+j*nr_spline_d, 
                D_Rat2_CUDA+j*nr_spline_d, nr_spline_d, r, f_cor);
    
    //if(f_cor<0)cout << "fcor  " << DiskDensf(r,0) << " " << Pot(r,0) << " " 
    //     << DiskDensfI(r,0,j) << " " << f_cor << endl;
    
    return DiskDensfICUDA(r,0,j, Radius_CUDA, A_Pot_CUDA)*f_cor;
}

__host__ __device__     
double RCircCUDA(double &am, thrust::device_ptr<double> Am_Tab_CUDA,
                 thrust::device_ptr<double> R_Tab_CUDA, 
                 thrust::device_ptr<double> R_Tab2_Zero_CUDA)
{
    double fam = fabs(am);
    
    if (fam > Am_Tab_CUDA[nr_d-1])
    {
        return R_Tab_CUDA[nr_d-1]*(fam/Am_Tab_CUDA[nr_d-1])*(fam/Am_Tab_CUDA[nr_d-1]);
    }
    else
    {
        double rcirc;
        SplintDCUDA(Am_Tab_CUDA, R_Tab_CUDA, R_Tab2_Zero_CUDA, nr_d, fam, rcirc);
        rcirc *= sqrt(fam);
        return rcirc;
    }
}

__host__ __device__     
double SigZ2CUDA(double r, int &i, thrust::device_ptr<double> Radius_CUDA, 
                 thrust::device_ptr<double> Rad_Spline_CUDA,
                 thrust::device_ptr<double> A_Pot_CUDA, 
                 thrust::device_ptr<double> FSZ_Rat_CUDA,
                 thrust::device_ptr<double> SZ_Rat2_CUDA)
{
    double psi_zh = PotCUDA(r, Z_Disk[i], Radius_CUDA, A_Pot_CUDA);
    double psi_0 = PotCUDA(r, 0.0, Radius_CUDA, A_Pot_CUDA);
    double true_sig_z2 = (psi_zh-psi_0)/log(0.419974);
    double f_cor;
    
    SplintDCUDA(Rad_Spline_CUDA+i*nr_spline_d, FSZ_Rat_CUDA+i*nr_spline_d, 
                SZ_Rat2_CUDA+i*nr_spline_d, nr_spline_d, r, f_cor);
    //cout << "sigz   " << FSZ_Rat[0] << " " << FSZ_Rat[1] << " " << SZ_Rat2[0] << " " << SZ_Rat2[1] << endl;
    return f_cor*true_sig_z2;
}

__host__ __device__     
double SigR2CUDA(double r, int &i)
{
    return Sigma_0[i]*Sigma_0[i]*exp(-r/R_Sigma[i]);
}

