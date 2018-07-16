//This implementation is NOT to be used for large l

__device__ double LegendreCUDA(const int l, const double x)
{
    if (l<0 || x<-1 || x>1)
    {
        cout << "Error! out of range legendre arguments" << endl;
        exit(1);
    }
    else if (l==0)
    {
        return 1;
    }
    else if (l==1)
    {
        return x;
    }
    else if (l==2)
    {
        return 0.5*(3*x*x-1);
    }
    else if (x = 1.0)
    {
        return 1;
    }
    else if (x==-1)
    {
        if (l%2==1)
        {
            return -1;
        }
        else if (l%2==0)
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
        
        for (int ell=2; ell<l+1; ++ell)
        {
            p_ell = (x*(2*ell-1)*p_ellm1 - (ell-1)*p_ellm2)/ell;
            p_ellm2 = p_ellm1;
            p_ellm1 = p_ell;
        }
        
        return p_ell;
    }
}
    
