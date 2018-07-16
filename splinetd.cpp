// Interpolation splines. From Numerical Recipes and so will be replaced.

#include "galaxy.h"

void SplintD(vector<double> &xa, vector<double> &ya, vector<double> &y2a, int n,
             double x, double &y)
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
        if (xa.at(k) > x)
            khi = k;
        else
            klo = k;
    }
        
    h = xa.at(khi)-xa.at(klo);
	
    if (h <= 0)
    {
        cerr << "bad xa input in splint. h = 0 " << xa.at(khi) << " " 
             << xa.at(klo) << " " << khi << " " << klo << endl;
        return;
    }
        
    double a = (xa.at(khi)-x)/h;
    double b = (x-xa.at(klo))/h;
      
    y = a*ya.at(klo) + b*ya.at(khi) + ((a*a*a-a)*y2a.at(klo) + (b*b*b-b)*y2a.at(khi))*h*h/6;
        
    return;
}

void SplineD(vector<double> &x, vector<double> &y, int n, double yp1, 
             double ypn, vector<double> &y2)
{
    vector<double> u(n);
    double qn, un;

    if (yp1 > 1e30)
    {
        y2.at(0) = 0;
        u.at(0) = 0;
    }
    else
    {
        y2.at(0)=-0.5;
        u.at(0) = (3/(x.at(1)-x.at(0)))*((y.at(1)-y.at(0))/(x.at(1)-x.at(0)) - yp1);
    }
    
    for (int i = 1; i < n-1; ++i)
    {
        double sig = (x.at(i)-x.at(i-1))/(x.at(i+1)-x.at(i-1));
        double p = sig*y2.at(i-1) + 2;
        y2.at(i) = (sig-1)/p;
        u.at(i) = (6*((y.at(i+1)-y.at(i))/(x.at(i+1)-x.at(i)) - (y.at(i)-y.at(i-1))/(x.at(i)-x.at(i-1))) /
               (x.at(i+1)-x.at(i-1)) - sig*u.at(i-1))/p;
        //cout << i << " "<<y.at(i+1)<<" "<<y.at(i)<<" "<<y.at(i-1)<<"    ";
    }
    
    if (ypn > 1e30)
    {
        qn = 0;
        un = 0;
    }
    else
    {
        qn = 0.5;
        un = (3/(x.at(n-1)-x.at(n-2))) * (ypn-(y.at(n-1)-y.at(n-2))/(x.at(n-1)-x.at(n-2)));
    }
    
    y2.at(n-1) = (un-qn*u.at(n-2))/(qn*y2.at(n-2)+1);

    for (int i = n-2; i > -1; --i)
    {
        //cout << y2.at(i) << "   " << u.at(i) << " ";
        y2.at(i) = y2.at(i)*y2.at(i+1) + u.at(i);
        //cout << u.at(i) << " ";
    }//cout << endl;
    
    return;
}
