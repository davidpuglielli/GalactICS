//Halo density profile functions

#include "galaxy.h"

void HaloDensityProfile(double &r, double *profile)
{
    double HaloProfileDens1 = HaloProfileDens(r);
    double dHaloProfileDens = HaloProfileDensPrime(r);
    double ddHaloProfileDens = HaloProfileDens2Prime(r);
    
    double trunc_fac = GetHaloTrunc(r);
    double dtrunc_fac = GetHaloTruncPrime(r);
    double ddtrunc_fac = GetHaloTrunc2Prime(r);
    
    profile[0] = HaloProfileDens1*trunc_fac;
    profile[1] = HaloProfileDens1*dtrunc_fac + dHaloProfileDens*trunc_fac;
    profile[2] = ddHaloProfileDens*trunc_fac + 2*dHaloProfileDens*dtrunc_fac + HaloProfileDens1*ddtrunc_fac;
}

double HaloDensity(double &r)
{
    double HaloProfileDens1 = HaloProfileDens(r);
    
    double trunc_fac = GetHaloTrunc(r);
    
    //cout << r << " " << HaloProfileDens1 << endl;
    
    return HaloProfileDens1*trunc_fac;
}

double HaloDensityPrime(double &r)
{
    double HaloProfileDens1 = HaloProfileDens(r);
    double dHaloProfileDens = HaloProfileDensPrime(r);
    
    double trunc_fac = GetHaloTrunc(r);
    double dtrunc_fac = GetHaloTruncPrime(r);
    
    return HaloProfileDens1*dtrunc_fac + dHaloProfileDens*trunc_fac;
}

double HaloDensity2Prime(double &r)
{
    double HaloProfileDens1 = HaloProfileDens(r);
    double dHaloProfileDens = HaloProfileDensPrime(r);
    double ddHaloProfileDens = HaloProfileDens2Prime(r);
    
    double trunc_fac = GetHaloTrunc(r);
    double dtrunc_fac = GetHaloTruncPrime(r);
    double ddtrunc_fac = GetHaloTrunc2Prime(r);
    
    //cout << "h2prime  " << HaloProfileDens << " " << dHaloProfileDens << " " << ddHaloProfileDens << " "
    //     << trunc_fac << " " << dtrunc_fac << " " << ddtrunc_fac << endl;
    //cout <<  ddHaloProfileDens*trunc_fac + 2*dHaloProfileDens*dtrunc_fac + HaloProfileDens*ddtrunc_fac << endl;
    
    return ddHaloProfileDens*trunc_fac + 2*dHaloProfileDens*dtrunc_fac + HaloProfileDens1*ddtrunc_fac;
}
    
double GetHaloTrunc(double &r)
{
    double truncfac;
    double erfarg = (r - G.c_halo) * oneoversqrt2 / G.drtrunc_halo;
    
    if (erfarg > 6) return truncfac = 0;
    else if (erfarg < -6) return truncfac = 1;
    else return truncfac = 0.5 * erfc(erfarg);
}

double GetHaloTruncPrime(double &r)
{
    double truncfac;
    double t = (r - G.c_halo) * oneoversqrt2 / G.drtrunc_halo;
    t *= t;
    
    if (t > 36) return truncfac = 0;
    else return truncfac = -exp(-t)/G.drtrunc_halo*oneoversqrt2pi;
}

double GetHaloTrunc2Prime(double &r)
{
    double truncfac;
    double t = (r - G.c_halo) * oneoversqrt2 / G.drtrunc_halo;
    double tt = t*t;
    
    if (tt > 36) return truncfac = 0;
    else return truncfac = t*exp(-tt)*oneoversqrtpi/G.drtrunc_halo/G.drtrunc_halo;
}

double HaloProfileDens(double &r)
{
//     if (r==0) r = 0.00001;
//     
//     double s = r/G.a_halo;
//     
//     double density = halo_const / pow(s,G.cusp) / pow(1+s,3-G.cusp);
//     
//     return density;
// 
//     double s = r/G.a_halo;
//     
//     //v_halo is actually rho_h here
//     double density = G.v_halo/(1+s*s);
//     
//     return density;
    
    if (!contraction_flag)
    {
        return RawHaloProfile(r);
    }
    else if (contraction_flag)
    {
        int rad_index = 0;//cout << r << endl;

        if (r <= Halo_AC_Radius[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius[nr_ac-1])
        {
            return RawHaloProfile(r);//*(1-baryon_frac);
        }
        else
        {
            while (Halo_AC_Radius[rad_index] < r)
            {
                ++rad_index;
            }
        }
        
        //cout << r << " " << rad_index << " " << endl;

        if (rad_index < 1)
        {
            cout << "Halo AC finds out of range indices. Exiting..." << endl;
            cout << "ihi " << rad_index << " " << r << " " << Halo_AC_Radius[1] << endl;
            exit(1);
        }
        else if (rad_index > nr_ac-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac-1;
        }

        double r1 = Halo_AC_Radius[rad_index-1];
        double r2 = Halo_AC_Radius[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << "Ac " << rad_index << " " << r << " " << t << " " << Halo_AC_Dens[rad_index]
        //     << "     " << tm1 << " " << Halo_AC_Dens[rad_index-1] << endl;
        
        return t*Halo_AC_Dens[rad_index] + tm1*Halo_AC_Dens[rad_index-1];     
    }       
}

double HaloProfileDensPrime(double &r)
{
//     if (r==0) r = 0.00001;
//     
//     double s = r/G.a_halo;
//     
//     double ddensity = -HaloProfileDens(r)*(3*s + G.cusp)/G.a_halo/s/(1 + s);
//     
//     return ddensity;
//     
//     double density = -HaloProfileDens(r)*2.0*r/(G.a_halo*G.a_halo+r*r);
//     
//     return density;
    
    //cout << "called" << endl;
    if (!contraction_flag)
    {
        return RawHaloProfilePrime(r);
    }
    else if (contraction_flag)
    {
        //cout << "called2" << endl;
        int rad_index = 0;

        if (r <= Halo_AC_Radius[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius[nr_ac-1])
        {
            return RawHaloProfilePrime(r);
        }
        else
        {
            while (Halo_AC_Radius[rad_index] < r)
            {
                ++rad_index;
            }
        }

        if (rad_index < 1)
        {
            cout << "Halo AC prime finds out of range indices. Exiting..." << endl;
            cout <<"ihi "<<rad_index<<" "<<r<< endl;
            exit(1);
        }
        else if (rad_index > nr_ac-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac-1;
        }

        double r1 = Halo_AC_Radius[rad_index-1];
        double r2 = Halo_AC_Radius[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        //cout << "Ac " << rad_index << " " << r << " " << t << " " << Halo_AC_Dens_D[rad_index]
        //     << "     " << tm1 << " " << Halo_AC_Dens_D[rad_index-1] << endl;
        
        return t*Halo_AC_Dens_D[rad_index] + tm1*Halo_AC_Dens_D[rad_index-1];
    }            
}

double HaloProfileDens2Prime(double &r)
{
//     if (r==0) r = 0.00001;
//     
//     double s = r/G.a_halo;
//     
//     double deno = s*s*G.a_halo*G.a_halo*(1+s)*(1+s);
//     double dddensity = HaloProfileDens(r)*(G.cusp*(G.cusp+1) + 8*G.cusp*s + 12*s*s) / deno;
//     
//     return dddensity;
    
//     double temp = G.a_halo*G.a_halo+r*r;
//     double density = (temp-4*r*r)/pow(temp,3);
//     density *= -2*G.v_halo*G.a_halo*G.a_halo;
//     
//     return density;
    
    if (!contraction_flag)
    {
        return RawHaloProfile2Prime(r);
    }
    else if (contraction_flag)
    {
        int rad_index = 0;

        if (r <= Halo_AC_Radius[1])
        {
            rad_index = 1;
        }
        else if (r > Halo_AC_Radius[nr_ac-1])
        {
            return RawHaloProfile2Prime(r);
        }
        else
        {
            while (Halo_AC_Radius[rad_index] < r)
            {
                ++rad_index;
            }
        }

        if (rad_index < 1)
        {
            cout << "Halo AC 2prime finds out of range indices. Exiting..." << endl;
            cout <<"ihi "<<rad_index<<" "<<r<< endl;
            exit(1);
        }
        else if (rad_index > nr_ac-1)
        {
            //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
            //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
            rad_index = nr_ac-1;
        }

        double r1 = Halo_AC_Radius[rad_index-1];
        double r2 = Halo_AC_Radius[rad_index];
        double t = (r-r1)/(r2-r1);
        double tm1 = 1-t;

        return t*Halo_AC_Dens_DD[rad_index] + tm1*Halo_AC_Dens_DD[rad_index-1];  
    }          
}

double RawHaloProfile(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/G.a_halo;
    
    double density = halo_const / pow(s,G.cusp) / pow(1+s,3-G.cusp);
    
    return density;
    
//     double s = r/G.a_halo;
//     double arg = -2.0*(pow(s, G.cusp)-1.0)/G.cusp;
//     
//     double density = G.v_halo*exp(arg);
//     
//     //cout << "          " << r << " " << density << endl;
//     
//     return density;
}

double RawHaloProfilePrime(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/G.a_halo;
    
    double ddensity = -HaloProfileDens(r)*(3*s + G.cusp)/G.a_halo/s/(1 + s);
    
    return ddensity;
    
//     double s = r/G.a_halo;
//     double arg = -2.0*(pow(s, G.cusp)-1.0)/G.cusp;
//     
//     double ddensity = -2*G.v_halo*exp(arg)*pow(r, G.cusp-1.0)/pow(G.a_halo, G.cusp);
//     
//     return ddensity;
}

double RawHaloProfile2Prime(double &r)
{
    if (r==0) r = 0.00001;
    
    double s = r/G.a_halo;
    
    double deno = s*s*G.a_halo*G.a_halo*(1+s)*(1+s);
    double dddensity = HaloProfileDens(r)*(G.cusp*(G.cusp+1) + 8*G.cusp*s + 12*s*s) / deno;
    
    return dddensity;
    
//     double s = r/G.a_halo;
//     double arg = -2.0*(pow(s, G.cusp)-1.0)/G.cusp;
//     double exparg = exp(arg);
//     double recip = 1.0/pow(G.a_halo, G.cusp);
//     double temp1 = (G.cusp-1.0)*pow(r, G.cusp-2.0);
//     double temp2 = 2.0*recip*pow(r, 2.0*G.cusp-2.0);
//     
//     double dddensity = -2*G.v_halo*exparg*recip*(temp1-temp2);
//     
//     return dddensity;
}

double GetHaloConst(void)
{
    return G.v_halo*G.v_halo*pow(2,1.0-G.cusp)*oneover4pi/G.a_halo/G.a_halo;
}
