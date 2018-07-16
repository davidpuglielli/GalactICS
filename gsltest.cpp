#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define ntheta 100

using namespace::std;

int main(void)
{
    for(int j=0; j<20000; ++j)
    {
        for (int l = 0; l < 12; ++l)
        {
            for (int i = 0; i < ntheta; ++i)
            {
                double costheta = double(i)/ntheta;
                gsl_sf_legendre_Pl(l, costheta);
                //pow(costheta,l);//cout<<costheta<<" "<<l<<endl;
            } 
        }
    }
    
    return 0;
}

