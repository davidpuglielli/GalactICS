#include "galaxy.h"

double GetTrunc(double &r, int &i)
{
    double truncfac;
    double erfarg = (r - G.Out_Disk[i]) * oneoversqrt2 / G.Dr_Trunc[i];
    
    if (erfarg > 4) return truncfac = 0;
    else if (erfarg < -4) return truncfac = 1;
    else return truncfac = 0.5 * erfc(erfarg);
}

void GetTruncPrime(double &r, int &i, double &truncfac, double &truncfacprime)
{
    double erfarg = (r - G.Out_Disk[i]) * oneoversqrt2 / G.Dr_Trunc[i];
    
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
        truncfacprime = -exp(-erfarg*erfarg) / G.Dr_Trunc[i] * oneoversqrt2pi;
    }
}

double GetTruncGas(double &r, int &i)
{
    double truncfac;
    double erfarg = (r - G.Out_GasDisk[i]) * oneoversqrt2 / G.Dr_Trunc_Gas[i];
    
    if (erfarg > 4) return truncfac = 0;
    else if (erfarg < -4) return truncfac = 1;
    else return truncfac = 0.5 * erfc(erfarg);
}

void GetTruncPrimeGas(double &r, int &i, double &truncfac, double &truncfacprime)
{
    double erfarg = (r - G.Out_GasDisk[i]) * oneoversqrt2 / G.Dr_Trunc_Gas[i];
    
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
        truncfacprime = -exp(-erfarg*erfarg) / G.Dr_Trunc_Gas[i] * oneoversqrt2pi;
    }
}

