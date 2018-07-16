//Get the distribution function from the integrals of motion (diskdf3), or from
//the velocity and position by calculating the integrals of motion (diskdf5)

#include "galaxy.h"

double DiskDF3ez(double &ep, double &am, double &ez, int &i)
{
    double psi00 = Pot(0, 0);
    double rc = RCirc(am);
    double freq_omega, freq_kappa;
    
    GetOmegaKappa(rc, freq_omega, freq_kappa);
    
    double v_circ = rc*freq_omega;
    
    double ec = -Pot(rc, 0) + 0.5*v_circ*v_circ;
    
    if (am < 0)
    {
        ec = -2*psi00-ec;
    }
    
    double sr2 = SigR2(rc,i), sz2 = SigZ2(rc,i);
    
    if (sz2 > 0)
    {
        double arg = (ec-ep)/sr2;
        if (arg>0)
            arg=0;
        double f_vert = FnaMidDen(rc,i)*exp(-ez/sz2)*oneoversqrt2pi/sqrt(sz2);
        f_vert = f_vert*freq_omega/freq_kappa/sr2*exp(arg)*oneoverpi;
        
        if(f_vert<0)
        {
            cout << "fvert < 0  " << freq_omega/freq_kappa/sr2 << " " 
                 << FnaMidDen(rc,i)<<endl;
            return 0;
        }
        
        return f_vert;
    }
    else
    {
        return 0;
    }
}
    
double DiskDF3intez(double &ep, double &am, double &ez, int &i)
{
    double psi00 = Pot(0, 0);
    double rc = RCirc(am);
    double freq_omega, freq_kappa;
    
    GetOmegaKappa(rc, freq_omega, freq_kappa);
    
    double v_circ = rc*freq_omega;
    
    double ec = -Pot(rc, 0) + 0.5*v_circ*v_circ;
    
    if (am < 0)
    {
        ec = -2*psi00-ec;
    }
    
    double sr2 = SigR2(rc,i), sz2 = SigZ2(rc,i);
    //if(rc>70)
        //cout << "df3   " << sz2 << " " << sr2 << " " << rc << " " << ec << "     " << freq_omega << " " 
        // << freq_kappa << " " << ep << " " << am << " " << ez << endl;
    //if(rc>70)
        //cout << "df3   " << sz2 << " " << freq_omega << " " << freq_kappa << endl;
    
    if (sz2 > 0)
    {
        double arg = (ec-ep)/sr2;
        if (arg>0)
            arg=0;
        double f_vert = FnaMidDen(rc,i)*exp(-ez/sz2);
                
        return freq_omega/freq_kappa/sqrt(sr2)*exp(arg)*f_vert*sqrt2*oneoversqrtpi;
    }
    else
    {
        return 0;
    }
}

double DiskDF5ez(double &vr, double &vt, double &vz, double &r, double &z, int &i)
{
    double psi_r0 = Pot(r, 0), psi_rz;
    
    if (z == 0)
    {
        psi_rz = psi_r0;
    }
    else
    {
        psi_rz = Pot(r, z);
    }
    
    double ep = 0.5*(vr*vr+vt*vt) - psi_r0;
    double am = r*vt;
    double ez = 0.5*vz*vz - psi_rz + psi_r0;
    
    double d = DiskDF3ez(ep, am, ez, i);
    
    //cout << "df   " << d << " " << ep << " " << am << " " << ez << endl;
    
    return d;
}

double DiskDF5intez(double &vt, double &r, double z, int &i)
{
    double psi_r0 = Pot(r, 0), psi_rz;
    
    if (z == 0)
    {
        psi_rz = psi_r0;
    }
    else
    {
        psi_rz = Pot(r, z);
    }
    
    double energy_plane = 0.5*(vt*vt) - psi_r0;
    double angular_momentum = r*vt;
    double energy_z = -psi_rz + psi_r0;
    
    double d = DiskDF3intez(energy_plane, angular_momentum, energy_z, i);
    
    //if(r>72)
        //cout << "df      " << psi_r0 << " " << psi_rz << " " << vt << " " << r << endl;
    //if(r>72)
        //cout << "dfint   " << d << " " << energy_plane << " " << angular_momentum << " " << energy_z << endl;
    
    return d;
}

 
double DiskDF3ezGas(double &ep, double &am, double &ez, int &i)
{
    double psi00 = Pot(0, 0);
    double rc = RCirc(am);
    double freq_omega, freq_kappa;
    
    GetOmegaKappa(rc, freq_omega, freq_kappa);
    
    double v_circ = rc*freq_omega;
    
    double ec = -Pot(rc, 0) + 0.5*v_circ*v_circ;
    
    if (am < 0)
    {
        ec = -2*psi00-ec;
    }
    
    double sr2 = SigR2Gas(rc,i), sz2 = SigZ2Gas(rc,i);
    
    if (sz2 > 0)
    {
        double f_vert = FnaMidDenGas(rc,i)*exp(-ez/sz2)*oneoversqrt2pi/sqrt(sz2);
        f_vert = f_vert*freq_omega/freq_kappa/sr2*exp((ec-ep)/sr2)*oneoverpi;
        
        return f_vert;
    }
    else
    {
        return 0;
    }
}
    
double DiskDF3intezGas(double &ep, double &am, double &ez, int &i)
{
    double psi00 = Pot(0, 0);
    double rc = RCirc(am);
    double freq_omega, freq_kappa;
    
    GetOmegaKappa(rc, freq_omega, freq_kappa);
    
    double v_circ = rc*freq_omega;
    
    double ec = -Pot(rc, 0) + 0.5*v_circ*v_circ;
    
    if (am < 0)
    {
        ec = -2*psi00-ec;
    }
    
    double sr2 = SigR2Gas(rc,i), sz2 = SigZ2Gas(rc,i);
    
    if (sz2 > 0)
    {
        double arg = (ec-ep)/sr2;
        if (arg>0)
            arg=0;
        double f_vert = FnaMidDenGas(rc,i)*exp(-ez/sz2);
        
        return freq_omega/freq_kappa*exp(arg)*f_vert/sqrt(sr2)*sqrt2*oneoversqrtpi;
    }
    else
    {
        return 0;
    }
}

double DiskDF5ezGas(double &vr, double &vt, double &vz, double &r, double &z, int &i)
{
    double psi_r0 = Pot(r, 0), psi_rz;
    
    if (z == 0)
    {
        psi_rz = psi_r0;
    }
    else
    {
        psi_rz = Pot(r, z);
    }
    
    double ep = 0.5*(vr*vr+vt*vt) - psi_r0;
    double am = r*vt;
    double ez = 0.5*vz*vz - psi_rz + psi_r0;
    
    //cout << "df   " << " " << ep << " " << am << " " << ez << endl;
    double d = DiskDF3ezGas(ep, am, ez, i);
    
    //cout << "df   " << d << " " << ep << " " << am << " " << ez << endl;
    
    return d;
}

double DiskDF5intezGas(double &vt, double &r, double z, int &i)
{
    double psi_r0 = Pot(r, 0), psi_rz;
    
    if (z == 0)
    {
        psi_rz = psi_r0;
    }
    else
    {
        psi_rz = Pot(r, z);
    }
    
    double energy_plane = 0.5*(vt*vt) - psi_r0;
    double angular_momentum = r*vt;
    double energy_z = -psi_rz + psi_r0;
    
    double d = DiskDF3intezGas(energy_plane, angular_momentum, energy_z, i);
    
    //if(r>72)cout << "df      " << psi_r0 << " " << psi_rz << " " << vt << " " << r << endl;
    //if(r>72)cout << "dfint   " << d << " " << energy_plane << " " << angular_momentum << " " << energy_z << endl;
    
    return d;
}

 
