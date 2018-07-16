//Get the halo density from the distribution function, using the energy as input

#include "galaxy.h"

double DensrPsi(double &rad, double &psi)
{
    int m = 33;//why?
    double v_min = -10;//actually log vmin

    if (psi < 0)
    {
        return 0;
    }
    
    double v_max = log10(sqrt(2*psi));
    double dlogv = (v_max - v_min) / (m-1);
    double logv, vel, energy, rho = 0;
    
    for (int i = 0; i < 4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        rho += 4*PI*dlogv*vel*vel*vel*DF(energy)*coef(i+1);
    }
    
    for (int i = 4; i < m-4; ++i)
    {
        logv = v_min + i*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        rho += 4*PI*dlogv*vel*vel*vel*DF(energy);
    }
    for (int j = 2; j < 5; ++j)
    {
        int jj = m-j+1;
        logv = v_min + jj*dlogv;
        vel = exp(logv);
        energy = psi - vel*vel/2;
        rho += 4*PI*dlogv*vel*vel*vel*DF(energy)*coef(j);
    }
    
    return rho;
}

double DF(double &energy)
{
    double dist_func, v02 = G.v_halo*G.v_halo, psi00 = Pot(0,0);
    double ev02 = energy/v02, psiv02 = 1 + psi00/v02;
    
    if (energy < 0 || ev02 > 1-psiv02)
    {
        return dist_func = 0;
    }
    else
    {
        double f1 = pow(ev02, 1.5);
        double f2 = pow(1-psiv02-ev02, -2.5);
        double ff3 = (1-psiv02-ev02) / (-log(ev02+psiv02));
        double f3 = pow(ff3, 2.715);
        double f4 = 0.362*ev02 - 0.5639*ev02*ev02 - 0.0859*ev02*ev02*ev02 - 
                    0.4912*pow(ev02,4);
        
        return dist_func = 0.09197*G.a_halo*f1*f2*f3*exp(f4);
    }
}

double coef(const int j)
{
    if (j == 1) return 17.0/48.0;
    if (j == 2) return 59.0/48.0;
    if (j == 3) return 43.0/48.0;
    if (j == 4) return 49.0/48.0;
    else 
    {
        cerr << "Bad coef input!" << endl;
        return -1;
    }
}
