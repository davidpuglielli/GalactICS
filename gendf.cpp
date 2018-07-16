//Get the distribution functions for the NFW halo and the Sersic bulge
//First the functions to generate the table, then the functions that calculate 
//the actual DF given the energy

#include "galaxy.h"

void GenTableE(void)
{
    double r = 0, r_min = 0.001, r_out = G.c_halo + 5*G.drtrunc_halo;

    psi_0 = GetTotalPsi(r);        cout << "    Psi_0    = " << psi_0 << endl;
    psi_crit = GetTotalPsi(r_out); cout << "    Psi_crit = " << psi_crit << endl;
    psi_d = GetTotalPsi(r_min);    cout << "    Psi_d    = " << psi_d << endl;
    
    for (int i = 0; i < n_psi; ++i)
    {
        double table_e = float(i)*log((psi_0-psi_crit)/(psi_0-psi_d))/float(n_psi-1);
        Table_E[i] = psi_0 - exp(table_e)*(psi_0-psi_d);
        
        //So that the halo edge is not artificially truncated by the else
        //statement in GenNFWDistFunc
        if (i == n_psi-1)
        {
            if (Table_E[i] != psi_crit)
            {
                Table_E[i] *= 1.0000000001;
            }
        }

        double start_rad = 1;
        
        psi_lower_initial = GetTotalPsi(start_rad);
        psi_upper_initial = GetTotalPsi(start_rad);
        log_rj2 = 1/log((psi_0 - psi_crit)/(psi_0 - psi_d));
                
        //cout << setw(12) << "Gentable " << psi_0 << " " << psi_crit << " "
        //     << psi_d << " " << Table_E[i] << " " << table_e << " " 
        //     << (psi_0-psi_crit)/(psi_0-psi_d) << " " << i << endl;
    }
}

void GenNFWDistFunc(void)
{
    double DF_NFW_last = -1000;
    
    ofstream DFFile;
    
    if (do_file_io)
    {
        DFFile.open("dfnfw.dat", ios::out);
    }
    
    for (int i = 0; i < n_psi; ++i)
    {
        double energy = Table_E[i];//cout<<"gennfwdf "<<i<<" "<<energy<<" "<<endl;
        DF_NFW[i] = GetNFWDistFunc(energy);//cout<<"DF   "<<i<<" "<<energy<<" "<<DF_NFW[i]<<endl;
        
        if (DF_NFW[i] > 0)
        {
            DF_NFW[i] = log(DF_NFW[i]);
            DF_NFW_last = DF_NFW[i];
        }
        else
        {
            DF_NFW[i] = DF_NFW_last;
            cout << "Warning: DF_NFW < 0 at bin " << i << ". Using previous "
                 << "value for the DF. This may\nyield the wrong mass and "
                 << "tidal radius." << endl;
        }
        
        if (do_file_io)
        {
            DFFile << setprecision(12) << Table_E[i] << "  " << DF_NFW[i] << endl;
        }
    }
}

double GetNFWDistFunc(double &energy)
{
    double t, t_max = sqrt(energy-psi_crit), dt = t_max/(n_int - 1);
    double psi = energy;
    double sum = 0, r_upper, r_lower;
    
    FindBrackets(psi, r_lower, r_upper);
    double tolerance = tolerance_factor;
    double r_psi = RootBisection(psi, r_upper, r_lower, tolerance);
    
    double d2rho_dpsi2 = GetNFWd2rhodpsi2(r_psi);
    sum += dt*d2rho_dpsi2;
    
    for (int i = 0; i < n_int-2; ++i)
    {
        t = (i + 1)*dt;
        psi = energy -t*t;
        FindBrackets(psi, r_lower, r_upper);
        tolerance = tolerance_factor;
        r_psi = RootBisection(psi, r_upper, r_lower, tolerance);
        d2rho_dpsi2 = GetNFWd2rhodpsi2(r_psi);
        
        sum += 2*dt*d2rho_dpsi2;
    }

    return sum*oneoversqrt8pi2;
}

void GenSersicDistFunc(void)
{
    double DF_Sersic_last = -1000;
    
    ofstream DFFile;
    
    if (do_file_io)
    {
        DFFile.open("dfsersic.dat", ios::out);
    }
    
    for (int i = 0; i < n_psi; ++i)
    {
        double energy = Table_E[i];//cout<<"gensersicdf "<<i<<" "<<energy<<" "<<endl;
        DF_Sersic[i] = GetSersicDistFunc(energy);
        
        if (DF_Sersic[i] > 0)
        {
            DF_Sersic[i] = log(DF_Sersic[i]);
            DF_Sersic_last = DF_Sersic[i];
        }
        else
        {
            DF_Sersic[i] = DF_Sersic_last;
            cout << "Warning: DF_Sersic < 0 at bin " << i << ". Using previous "
                 << "value for the DF. This may\nyield the wrong mass and "
                 << "tidal radius." << endl;
        }
        
        if (do_file_io)
        {
            DFFile << setprecision(12) << Table_E[i] << "  " << DF_Sersic[i] << endl;
        }
    }
    
}

double GetSersicDistFunc(double &energy)
{
    double t, t_max = sqrt(energy-psi_crit), dt = t_max/(n_int - 1);
    double DF, psi = energy;
    double sum = 0, r_upper, r_lower;
    
    FindBrackets(psi, r_lower, r_upper);
    double tolerance = tolerance_factor;
    double r_psi = RootBisection(psi, r_upper, r_lower, tolerance);
    
    double d2rho_dpsi2 = GetSersicd2rhodpsi2(r_psi);
    sum += dt*d2rho_dpsi2;
    
    for (int i = 0; i < n_int-2; ++i)
    {
        t = (i + 1)*dt;
        psi = energy - t*t;
        FindBrackets(psi, r_lower, r_upper);
        tolerance = tolerance_factor;
        r_psi = RootBisection(psi, r_upper, r_lower, tolerance);
        d2rho_dpsi2 = GetSersicd2rhodpsi2(r_psi);
        
        sum += 2*dt*d2rho_dpsi2;
    }
    
    return DF = sum*oneoversqrt8pi2;
}

//what about the cusp?
//Wasteful - we calculate the force twice when we should only need to once. Fix!
double GetNFWd2rhodpsi2(double &r)
{
    double density = HaloDensity(r);//cout<<"rhooo   "<<r<<" "<<density<<endl;
    double ddensity = HaloDensityPrime(r);//cout<<"rhoop   "<<r<<" "<<ddensity<<endl;
    double dddensity = HaloDensity2Prime(r);//cout<<"rhopp   "<<r<<" "<<dddensity<<endl;
    double force = HaloForce(r);//cout<<"force   "<<r<<" "<<force<<endl;
    
    if (bulge_flag == 1)
    {
        force += SersicForce(r);
        density += SersicDens(r);
    }
    
    if (disk_flag == 1)
    {
        force += DiskForce(r);
        density += DiskDensity(r);
    }
    
    double t1 = fourpi*density*ddensity/force;
    double t2 = 2*ddensity/r;
    return (t1 + t2 + dddensity)/(force*force);
}

double GetSersicd2rhodpsi2(double &r)
{
    double density = SersicDens(r);
    double ddensity = SersicDensPrime(r);
    double dddensity = SersicDens2Prime(r);
    double force = SersicForce(r);
    
    if (halo_flag == 1)
    {
        force += HaloForce(r);
        density += HaloDensity(r);
    }
    
    if (disk_flag == 1)
    {
        force += DiskForce(r);
        density += DiskDensity(r);
    }
    
    double t1 = fourpi*density*ddensity/force;
    double t2 = 2*ddensity/r;
    
    return (t1 + t2 + dddensity)/(force*force);
}
    
void FindBrackets(double &psi, double &r_lower, double &r_upper)
{
//    r_lower=Radius[0];
//    r_upper=r_max;//cout<<psi<<endl;
     r_lower = r_upper = 1;//sqrt(r_max*dr);
     
     double psi_lower = psi_lower_initial;//GetTotalPsi(r_lower);
     double psi_upper = psi_upper_initial;//GetTotalPsi(r_upper);
     
     //cout<<"findbrackets "<<psi << " " <<psi_upper<<" "<<psi_lower<<endl;
     
     while(psi_lower <= psi)
     {
         r_lower*=0.1;
         psi_lower=GetTotalPsi(r_lower);
         //cout<<" findbrackets1 "<<psi << " " <<r_lower<<" " <<psi_upper<<" "<<psi_lower<<endl;
     }
     
     while (psi_upper >= psi)
     {
         r_upper=r_upper*10;//min(r_upper*10,r_max);
         psi_upper=GetTotalPsi(r_upper);
         //cout<<" findbrackets2 "<<psi << " "<<r_upper<<" " <<psi_upper<<" "<<psi_lower<<endl;
     }
     
     //cout<<"findbrackets   "<<r_upper<<" "<<r_lower<<endl;
}

double RootBisection(double &psi, double &r_upper, double &r_lower, double tolerance)
{
    double psi_upper = GetTotalPsi(r_upper)-psi;
    double psi_lower = GetTotalPsi(r_lower)-psi;
    double r_midpoint, psi_midpoint, r_gap;
    int iters = 0;
    
    //cout << "Bisection  " << psi << " " << psi_upper+psi << " " << psi_lower+psi << endl;
    
    if (psi_upper*psi_lower > 0)
    {
        cerr << "Bisection endpoints do not bracket root. Exiting..." << endl;
        cerr << psi_upper << " " << psi_lower << " " << r_lower << " " 
             << r_upper << endl;
        exit(1);
    }
    
    do
    {
        ++iters;
        //cout<<"doing stuff "<<max_iter<<" "<<r_gap<<" "<<iter<<" "<<tolerance<<endl;
        
        r_midpoint = (r_upper+r_lower)/2;
        
        psi_midpoint = GetTotalPsi(r_midpoint)-psi;
        
        if (psi_lower*psi_midpoint < 0)
        {
            r_upper = r_midpoint;
            psi_upper = psi_midpoint;
        }
        else if (psi_upper*psi_midpoint < 0)
        {
            r_lower = r_midpoint;
            psi_lower = psi_midpoint;
        }
        else if (psi_midpoint==0)
        {
            return r_midpoint;
            cout << "\nFinished root bisection1\n" << endl;
        }
        else
        {
            cerr << "Root bisection failed! exiting..." << endl;
            exit(1);
        }
        
        r_gap = r_upper-r_lower;
    }
    while (iters < max_iter && r_gap > tolerance);
    
    //cout << "Finished root bisection " << iters << endl;
    
    return r_midpoint;
}
