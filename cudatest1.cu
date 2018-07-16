#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <iostream>

using namespace::std;

__device__ double varb = -10.3;
//double varb=2;

struct galstruct
{
    double i;
    double j;
    galstruct(double ii, double jj){i=ii;j=jj;}
};

__host__ __device__ double Pot(double s, double z, thrust::device_ptr<float> denspsihalo)
{
    //printf("%f\n", *(denspsihalo+12));
    return s*z*(denspsihalo[12]);
}

__host__ __device__ double Pot(double s, double z)
{
    //printf("%f\n", *(denspsihalo+12));
    return s*z;
}

__host__ __device__ double PotTrunc(galstruct G, double s, double z)
{
    return G.i*s*z;
}

struct saxpy_functor
{
    const galstruct a;
    //const float b;
    const thrust::device_ptr<float> Vec;
    const thrust::device_ptr<float> Vec2;
    //const vector<float> Vec;
    //const vector<float> Vec2;
    
    saxpy_functor(galstruct _a, thrust::device_ptr<float> _Vec, 
                  thrust::device_ptr<float> _Vec2):a(_a),Vec(_Vec),Vec2(_Vec2){}
    //saxpy_functor(galstruct _a, float _b):a(_a), b(_b){}
    //saxpy_functor(){}
    //saxpy_functor(galstruct _a, vector<float> _Vec, 
    //              vector<float> _Vec2):a(_a),Vec(_Vec),Vec2(_Vec2){}
            
    saxpy_functor(galstruct _a):a(_a){}
            
    __host__ __device__
        float operator()(const float& x, const float& y) const
        {
            double psi = 1.4;
            double r=x,costheta=y;
            double z=r*costheta, rad=r*sqrt(1-costheta*costheta);
            //if (y>a.i) return 2.1;
            
            double rt = (x+y)/2.0;
            //if (i > 999) i = 999;
            
            double trunc_fac=PotTrunc(a,x,y);
            
            double con;
            
            if (z==0)
            {
                con = 1;
            }
            else
            {
                //KK: make disk vertically isothermal; normalize density so that
                //at height 3*zdisk the vertical density drop is that of a sech^2.
                double psizh = Pot(2.0, 4, Vec2);//Pot(r, 3.0*G.Z_Disk[i]);
                double psi00 = 2.0;//Pot(r, 0);
                double dpsizh = psi00 - psizh;
                double dpsi = psi00 - psi;
                double coeff = 0;

                if (r==0 && z==0)
                {
                    dpsi = 0;
                }

                if (fabs(dpsizh) > 0)
                {
                    coeff = dpsi / dpsizh;
                }

                if (coeff > 16 || coeff < 0)
                {
                    con = psizh;
                }
                else
                {
                    con = pow(0.009866, coeff);
                }
            }
            
            //return disk_const*s_profile*con*trunc_fac;
            return a.i*x*y*(Vec[rt])*con;
        }
};

int main(void)
{
    //int major = THRUST_MAJOR_VERSION;
    //int minor = THRUST_MINOR_VERSION;
    
    //std::cout << "Thrust v" << major << "." << minor << std::endl;
    varb = -10.3;
    thrust::device_vector<float> D(1000);
    thrust::device_vector<float> E(1000, 3.14159);
    thrust::device_vector<float> F(1000);
    thrust::device_vector<float> G(1000, 2.71828);
    
    //float *rawptr;
    //cudaMalloc((void **)&rawptr, 1000*sizeof(float));
    //thrust::device_ptr<float> dev_ptr(rawptr);
    vector<double> Dens_Psi_Halo(1000,60);
    
    thrust::device_ptr<float> dev_ptr = thrust::device_malloc<float>(1000);
    thrust::device_ptr<float> dev_ptr2 = thrust::device_malloc<float>(1000);
    thrust::device_ptr<float> denspsihalo = thrust::device_malloc<float>(1000);
    
    thrust::sequence(D.begin(),D.end());
    thrust::sequence(dev_ptr, dev_ptr+1000);
    thrust::copy(Dens_Psi_Halo.begin(), Dens_Psi_Halo.end(), denspsihalo);    
    
/***    thrust::device_vector<float> dev_ptr(1000);
    thrust::device_vector<float> denspsihalo(1000);
    
    thrust::sequence(D.begin(),D.end());
    thrust::sequence(dev_ptr.begin(), dev_ptr.end());
    thrust::copy(Dens_Psi_Halo.begin(), Dens_Psi_Halo.end(), denspsihalo.begin());
***/    
        const galstruct A(10.2,2);
    //A.i = 0.2;
    //A.j = 2;
    
    //thrust::transform(D.begin(),D.end(),E.begin(),F.begin(),thrust::multiplies<float>());
    
    thrust::transform(D.begin(),D.end(),E.begin(),dev_ptr2,saxpy_functor(A,dev_ptr,denspsihalo));
    
    for (int i = 0; i < D.size(); ++i)
    {
        std::cout << i << " " << D[i] << " " << E[i] << " " << F[i] << " " << dev_ptr[i] << std::endl;
    }
    
    std::cout << varb << " " << PotTrunc(A, 3.2, 3) << std::endl;
    
    return 0;
}
