// Get the angular and epicycle frequencies as functions of r (GetEpicycleFrequencies)
// and at a radius r (GetOmegaKappa) by spline fitting those functions

#include "galaxy.h"

void GetOmegaKappa(double &r, double &freq_omega, double &freq_kappa)
{
    // KK: a anrd b are coefficients of potential -a/r -b/r^3 +c 
    // which give omeka and kappa as observed
    // at the outer radius of the model---these are used to extrapolate if needed.
    // in this potential om^2=a/r^3+3b/r^5, kap^2=a/r^3-3b/r^5 
    double rad = Radius[nr-1];
    double a = rad*rad*rad/2*(Omega[nr-1]*Omega[nr-1]+A_K[nr-1]*A_K[nr-1]);
    double b = pow(rad, 5)/6*(Omega[nr-1]*Omega[nr-1]-A_K[nr-1]*A_K[nr-1]);
    //double c = a/rad/rad/rad+b/pow(rad, 5)+Pot_Major_Tot[nr-1];
    
    //cout << "Getting OmeKap " << r << endl;
    
    if (r > Radius[nr-1])
    {
        double r3 = a/rad/rad/rad;
        double r5 = 3*b/pow(rad, 5);
        freq_omega = sqrt(r3+r5);
        freq_kappa = sqrt(r3-r5);
    }
    else
    {
        SplintD(Radius, Omega, Omega2, nr, r, freq_omega);
        SplintD(Radius, A_K, A_K2, nr, r, freq_kappa);
    }
    //cout << "Getting OmeKap " << r << endl;
}

void GetEpicycleFrequencies(void)
{
    // KK: read omega, potential anrd potential'' from the frequencies file,
    // spline-fit pot(r), and then use the 2nd derivatives of the
    // potential in the table of kappa. these can then in turn be
    // spline-fitted. should arrange for the boundary condition at r=0 to
    // be even for pot, omega and kappa.
    
    // DP: no extra file IO introduced here unless requested.
    
    //vector <double> Omega(nr), Omega2(nr), A_K(nr), A_K2(nr);
    ofstream omekap;
    
    cout << "Getting Kap" << endl;
    
    for (int i = 1; i < nr; ++i)
    {
        Omega[i] = VC_Tot[i]/Radius[i];
        //cout << i << " " << Omega[i] << endl;
    }
    
    // KK: extrapolate omega (linearly) and potential (quadratically) to zero
    Omega[0] = 2*Omega[1] - Omega[3];
    Pot_Major_Tot[0] = (4*Pot_Major_Tot[1]-Pot_Major_Tot[3])/3;
    Psi2[0] = (4*Psi2[1]-Psi2[3])/3;
    
    SplineD(Radius, Omega, nr, 1e32, 1e32, Omega2);
    
    // KK: calculate epicycle frequencies
    A_K[0] = 2*Omega[0];
    
    for (int i = 1; i < nr; ++i)
    {
        double temp = Psi2[i]+3*Omega[i]*Omega[i];
        A_K[i] = sqrt(max(0.0,temp));
        //cout << i << " " << A_K[i] << endl;
    }
    
    SplineD(Radius, A_K, nr, 0, 1e32, A_K2);
    
//     for (int i = 0; i < nr; ++i)
//     {
//         cout << Radius[i] << " " << Omega[i] << " " << Omega2[i] << " " << A_K[i] 
//              << " " << A_K2[i] << " " << Pot_Major_Tot[i] << " " << Psi2[i] << endl;
//     }
    
	if (do_file_io)
	{
        // KK: a anrd b are coefficients of potential -a/r -b/r^3 +c 
        // which give omeka and kappa as observed
        // at the outer radius of the model---these are used to extrapolate if needed.
        // in this potential om^2=a/r^3+3b/r^5, kap^2=a/r^3-3b/r^5 
        double rad = Radius[nr-1];
        double a = rad*rad*rad/2*(Omega[nr-1]*Omega[nr-1] + A_K[nr-1]*A_K[nr-1]);
        double b = pow(rad, 5)/6*(Omega[nr-1]*Omega[nr-1] - A_K[nr-1]*A_K[nr-1]);
        double c = a/rad/rad/rad + b/pow(rad, 5) + Pot_Major_Tot[nr-1];
    
	    omekap.open("omekap.dat");
		
		for (int i = 0; i < nr; ++i)
		{
		    omekap << Radius[i] << " " << Omega[i] << " " << Omega2[i] << " "
			       << A_K[i] << " " << A_K2[i] << " " << Pot_Major_Tot[i] << " "
			       << Psi2[i] << endl;
		}
	
        for (int i = 0; i < 10; ++i)
        {
            double radx = Radius[nr-1]*(1+0.1*i);
            double r3 = a/radx/radx/radx;
            double r5 = 3*b/pow(radx, 5);
            double freq_omega = sqrt(r3+r5);
            double freq_kappa = sqrt(r3-r5);
         
            omekap << radx << " " << freq_omega << " " << 0 << " " << freq_kappa << " "
                   << 0 << " " << c-(r3+r5/3)*radx*radx << " " << -2*r3-4*r5 << endl;
        }
    }
}

double FnaMidDenGas(double &r, int &j)
{
    double f_cor;//, dens_const = GetDensConst(r, j);
    
    SplintD(Rad_Spline_Gas[j], FD_Rat_Gas[j], D_Rat2_Gas[j], nr_spline, r, f_cor);
    
    //if(f_cor<0)cout << "fcor  " << DiskDensf(r,0) << " " << Pot(r,0) << " " 
    //     << DiskDensfI(r,0,j) << " " << f_cor << endl;
    
    return GasDiskDensfI2(r,0,j)*f_cor;
}

double FnaMidDen(double &r, int &j)
{
    double f_cor;
    
    SplintD(Rad_Spline[j], FD_Rat[j], D_Rat2[j], nr_spline, r, f_cor);
    
//     if(f_cor<0)
//     {
//         cout << "fcor  " << DiskDensf(r,0) << " " << Pot(r,0) << " " 
//              << DiskDensfI(r,0,j) << " " << f_cor << endl;
//         for (int i = 0; i < nr_spline; ++i)
//         cout << Rad_Spline[j][i]<<" "<< FD_Rat[j][i]<<" "<< D_Rat2[j][i]<<endl;
//     }
    
    return DiskDensfI(r,0,j)*f_cor;
}

void GetRCirc(void)
{
    // KK: make a spline fit to [Rcirc/sqrt(Lcirc)] vs. Lcirc.
    // this is ~ constant at the origin, and goes like Lcirc**1.5 
    // at large radii (where the potential is Kepler).
    
    //cout << "Getting RCirc " << endl;
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        double om, t;
        GetOmegaKappa(rad,om,t);
        Am_Tab[i] = om*rad*rad;
        R_Tab[i] = 1/sqrt(om);
        R_Tab2_Zero[i] = 0;
        //cout << "rc      " << i << " " << rad << " " << om << " " << Am_Tab[i] << " " << R_Tab[i] <<  endl;
    }
    
    double slope_inf = 1.5*R_Tab[nr-1]/Am_Tab[nr-1];
    
    SplineD(Am_Tab, R_Tab, nr, 1e32, slope_inf, R_Tab2);
    
    //cout << "Getting RCirc " << endl;
    //for (int i = 0; i < nr; ++i)
    //{
        //cout << "rcirc   " << i << " " << Am_Tab[i] << " " << R_Tab[i] << " " << R_Tab2[i] << endl;
    //}
}

double RCirc(double &am)
{
    double fam = fabs(am);
    
    if (fam > Am_Tab[nr-1])
    {
        return R_Tab[nr-1]*(fam/Am_Tab[nr-1])*(fam/Am_Tab[nr-1]);
    }
    else
    {
        double rcirc;
        SplintD(Am_Tab, R_Tab, R_Tab2_Zero, nr, fam, rcirc);
        rcirc *= sqrt(fam);
        return rcirc;
    }
}
