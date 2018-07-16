#include "galaxy.h"

void DiskDF(void)
{
    double freq_omega, freq_kappa;
    
    //Initialise spline fitting functions    
    GetEpicycleFrequencies();
    GetRCirc();
    
    ofstream cordbhout;
    
    if (do_file_io)
    {
        cordbhout.open("cordbh.dat");
    }
    
    for (int i = 0; i < disk; ++i)
    {
        double dr_spline = Disk_Edge[i]/nr_spline;
    
        for (int j = 0; j < nr_spline+1; ++j)
        {
            Rad_Spline[i][j] = j*dr_spline;
            //cout << Rad_Spline[j] << endl;
        }
    
        double rfid = 2.5*G.R_Disk[i];
        GetOmegaKappa(rfid, freq_omega, freq_kappa);
        
        cout << "omega   " << freq_omega << " " << freq_kappa << endl;
        
        double sigma_r = sqrt(SigR2(rfid,i));
        double sigden = DiskDensfI(rfid, 0, i)*2*G.Z_Disk[i];
        
        double sigma_r_crit = 3.36*sigden/freq_kappa;
        double Q_toomre = sigma_r/sigma_r_crit;
        
        cout << "Sigma_r_crit = " << sigma_r_crit << " " << sigden << " " 
             << freq_kappa << " " << DiskDensfI(rfid, 0, i) << " " 
             << freq_omega << " at R = 2.5R_d" << endl;
        cout << "Toomre Q = " << Q_toomre << " at R = 2.5R_d" << endl;
        
        // KK: first call initializes the spline fitting functions by the way
        // DP: That's now done in a separate function, GetEpicycleFrequencies().
        GetOmegaKappa(G.Out_Disk[i], freq_omega, freq_kappa);
        
        // DP: I still have no idea what this message means
        if (sqrt(SigR2(G.Out_Disk[i],i)) > freq_kappa*G.Dr_Trunc[i])
        {
            cout << " WARNING : You have asked for an outer velocity dispersion"
                 << " which is larger\n           than that used in the disc"
                 << " truncation. Strange things may happen." << endl;
            cout << sqrt(SigR2(G.Out_Disk[i],i))/(G.Dr_Trunc[i]*freq_kappa) << " " 
                 << freq_kappa << endl;
        }
        
        for (int j = 0; j < nr_spline+1; ++j)
        {
            FD_Rat[i][j] = 1;
            FSZ_Rat[i][j] = 1;
        }
        
        SplineD(Rad_Spline[i], FD_Rat[i], nr_spline+1, 1e32, 1e32, D_Rat2[i]);
        SplineD(Rad_Spline[i], FSZ_Rat[i], nr_spline+1, 1e32, 1e32, SZ_Rat2[i]);
        
        //for (int j = 0; j < nr_spline; ++j)
        //{
        //    cout << "Spline   " << Rad_Spline[j] << " " << FD_Rat[j] << " " 
        //         << D_Rat2[i] << " " << FSZ_Rat[j] << " " << SZ_Rat2[i] << endl;
        //}
        
        // KK: Make niter iterations on the correction functions for 
        // midplane density and vertical velocity dispersion.
        // DP: Integration is compared to model disk structure from the density
        // functions. The correction functions measure the deviation from the
        // target disk density given by DiskDensfI().
        for (int k = 0; k < n_iter; ++k)
        {
            for (int j = 0; j < nr_spline+1; ++j)
            {
                double rad = Rad_Spline[i][j];
                
                GetOmegaKappa(rad, freq_omega, freq_kappa);
                //cout << "freqs   " << k << " " << j << " " << freq_omega << " " 
                //     << freq_kappa << endl;
                
                double sig_r = sqrt(SigR2(rad,i));
                
                double v_circ = rad*freq_omega;
                
                double d0 = 0, dz2 = 0;
                double dvr = 0.1*sig_r;
        
                // KK: integrate analytically over vr and vz, numerically over vphi
                for (int ivt = 1; ivt < 102; ++ivt)
                {
                    double vt = v_circ+(ivt-51)*dvr;
                    double df = DiskDF5intez(vt, rad, 0, i);
                    double c = 0.33333333333333*((ivt%2)*2+2);
                    d0 += c*dvr*df;
                    double df1 = DiskDF5intez(vt, rad, G.Z_Disk[i], i);
                    dz2 += c*dvr*df1;
                }
                
                D_Rat[i][j] = d0/DiskDensfI(rad, 0, i);
                DZ2_Rat[i][j] = dz2/DiskDensfI(rad, G.Z_Disk[i], i);
                FZ_Rat[i][j] = log(DiskDensfI(rad, 0, i)/DiskDensfI(rad, G.Z_Disk[i], i))/log(d0/dz2);
            }
            
            //Output a table of residuals to the console
             cout << endl;
             
             for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
             {
                 cout << setprecision(5) << left << Rad_Spline[i][j] << " ";
             }
             
             cout << endl;
             cout << "0.000 ";
             
             for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
             {
                 cout << setprecision(5) << left << D_Rat[i][j] << " ";
             }
             
             cout << endl;
             cout << G.Z_Disk[i] << " ";
             
             for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
             {
                 cout << setprecision(5) << left << DZ2_Rat[i][j] << " ";
             }
             
             cout << endl;
            
            for (int j = 0; j < nr_spline+1; ++j)
            {
                FD_Rat[i][j] = FD_Rat[i][j]/D_Rat[i][j];
                FSZ_Rat[i][j] = FSZ_Rat[i][j]/FZ_Rat[i][j];
            }
            
            SplineD(Rad_Spline[i], FD_Rat[i], nr_spline+1, 1e32, 1e32, D_Rat2[i]);
            SplineD(Rad_Spline[i], FSZ_Rat[i], nr_spline+1, 1e32, 1e32, SZ_Rat2[i]);
        }
        
        // KK: Finally, write out a more detailed table of residuals.
        for (int j = 0; j < nr_spline+1; ++j)
        {
            double rad = Rad_Spline[i][j];

            GetOmegaKappa(rad, freq_omega, freq_kappa);

            double sig_r = sqrt(SigR2(rad,i));

            double v_circ = rad*freq_omega;

            double d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0;
            double dvr = 0.1*sig_r;

            for (int ivt = 1; ivt < 102; ++ivt)
            {
                double vt = v_circ+(ivt-51)*dvr;
                double df = DiskDF5intez(vt, rad, 0, i);
                double c = 0.33333333333333*((ivt%2)*2+2);
                d0 += c*dvr*df;
                df = DiskDF5intez(vt, rad, 0.5*G.Z_Disk[i], i);
                d1 += c*dvr*df;
                df = DiskDF5intez(vt, rad, G.Z_Disk[i], i);
                d2 += c*dvr*df;
                df = DiskDF5intez(vt, rad, 1.5*G.Z_Disk[i], i);
                d3 += c*dvr*df;
                df = DiskDF5intez(vt, rad, 2*G.Z_Disk[i], i);
                d4 += c*dvr*df;
            }
            
            dirat[0][j] = d0/DiskDensfI(rad, 0, i);
            dirat[1][j] = d1/DiskDensfI(rad, 0.5*G.Z_Disk[i], i);
            dirat[2][j] = d2/DiskDensfI(rad, G.Z_Disk[i], i);
            dirat[3][j] = d3/DiskDensfI(rad, 1.5*G.Z_Disk[i], i);
            dirat[4][j] = d4/DiskDensfI(rad, 2*G.Z_Disk[i], i);
        }
        
        //Output a table of residuals to the console
        cout << endl;
        
        for (int j = 0; j < nr_spline; j+=nr_spline/10)
        {
            cout << setprecision(5) << left << Rad_Spline[i][j] << " ";   
        }
        
        cout << endl;
        for (int k = 0; k < 5; ++k)
        {
            cout << 0.5*k*G.Z_Disk[i] << " ";
        
            for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
            {
                cout << setprecision(5) << left << dirat[k][j] << " ";
            }
        
            cout << endl;
        }
        
        cout << endl;
            
        //write to cordbh.dat for each disk component if requested
        if (do_file_io)
        {
            cordbhout << G.Sigma_0[i] << " " << G.R_Sigma[i] << " " << nr_spline << endl;
            
            for (int j = 0; j < nr_spline+1; ++j)
            {
                cordbhout << Rad_Spline[i][j] << " " << D_Rat[i][j] << " " 
                          << D_Rat2[i][j] << " " << FD_Rat[i][j] << " " 
                          << FZ_Rat[i][j] << " " << FSZ_Rat[i][j] << endl;
            }
        }
        else
        {
            //for (int j = 0; j < nr_spline+1; ++j)
            //{
            //    cout << setw(12) << Rad_Spline[i][j] << " " << setw(12) 
            //         << D_Rat[i][j] << " " << setw(12) << D_Rat2[i][j] << " "
            //         << setw(12) << FD_Rat[i][j] << " " << setw(12) 
            //         << FZ_Rat[i][j] << " " << setw(12) << FSZ_Rat[i][j] << endl;
            //}
        }
        
        //Now get the Toomre stability parameters Q and X
        double q_avg = 0, x_avg = 0, q_scale, x_scale;
        int points = 0;
        ofstream toomreq("toomre.dat"), toomrex("toomres_other_parameter.dat");
        
        for (int j = 0; j < nr_spline; ++j)
        {
            double rad = Rad_Spline[i][j];
            GetOmegaKappa(rad, freq_omega, freq_kappa);
            double sig_r = sqrt(SigR2(rad,i));
            double sig_den = DiskDensfI(rad,0,i)*2*G.Z_Disk[i];
            
            double mode = 2;
            double q = sig_r*freq_kappa/3.36/sig_den;
            double x = freq_kappa*freq_kappa*rad*oneover2pi/mode/sig_den;
            
            if (rad < 2*G.R_Disk[i])
            {
                q_scale = q;
                x_scale = x;
            }
            
            if (rad > 0.5*G.R_Disk[i] && rad < min(3*G.R_Disk[i], Disk_Edge[i]))
            {
                q_avg += q;
                x_avg += x;
                ++points;
            }
            
            toomreq << rad << " " << q << endl;
            toomrex << rad << " " << x << endl;
        }
            
        Q_scale.push_back(q_scale);
        X_scale.push_back(x_scale);
        Q_avg.push_back(q_avg/points);
        X_avg.push_back(x_avg/points);
            
        cout << "Q at 2 scale lengths = " << Q_scale[i] << endl;
        cout << "X at 2 scale lengths = " << X_scale[i] << endl;
        cout << "Q average = " << Q_avg[i] << endl;
        cout << "X average = " << X_avg[i] << endl;
        cout << endl;
        
    }

    //To calculate the gas disk DF, we start by assuming an ideal gas and use that
    //to calculate the temperature at any point. Then we calculate the dispersion
    //from the temperature.
    
    for (int i = 0; i < gas_disk; ++i)
    {
        cout << "Now getting gas disk DF" << endl;
        double dr_spline = GasDisk_Edge[i]/nr_spline;
    
        for (int j = 0; j < nr_spline+1; ++j)
        {
            Rad_Spline_Gas[i][j] = j*dr_spline;
            //cout << Rad_Spline[j] << endl;
        }
    
        double rfid = 2.5*G.R_GasDisk[i];
        GetOmegaKappa(rfid, freq_omega, freq_kappa);
        
        cout << "omega   " << freq_omega << " " << freq_kappa << endl;
        
        double sigma_r = sqrt(SigR2Gas(rfid,i));
        double sigden = GasDiskDensfI2(rfid, 0, i)*2*G.Z_GasDisk[i];
        
        double sigma_r_crit = PI*sigden/freq_kappa;
        double Q_toomre = sigma_r/sigma_r_crit;
        
        cout << "Sigma_r_crit = " << sigma_r_crit << " " << sigden << " " 
             << freq_kappa << " " << GasDiskDensfI2(rfid, 0, i) << " " 
             << freq_omega << " at R = 2.5R_d" << endl;
        cout << "Toomre Q = " << Q_toomre << " at R = 2.5R_d" << endl;
        
        // KK: first call initializes the spline fitting functions by the way
        // DP: That's now done in a separate function, GetEpicycleFrequencies().
        GetOmegaKappa(G.Out_GasDisk[i], freq_omega, freq_kappa);
        
        cout << "Now getting gas disk DF" << endl;
        // DP: I still have no idea what this message means
        if (sqrt(SigR2Gas(G.Out_GasDisk[i],i)) > freq_kappa*G.Dr_Trunc_Gas[i])
        {
            cout << " WARNING : You have asked for an outer velocity dispersion"
                 << " which is larger\n           than that used in the disc"
                 << " truncation. Strange things may happen." << endl;
            cout << sqrt(SigR2Gas(G.Out_GasDisk[i],i))/(G.Dr_Trunc_Gas[i]*freq_kappa) << " " 
                 << freq_kappa << endl;
        }
        
        cout << "Now getting gas disk DF" << endl;
        for (int j = 0; j < nr_spline+1; ++j)
        {
            FD_Rat_Gas[i][j] = 1;
            FSZ_Rat_Gas[i][j] = 1;
        }
        
        SplineD(Rad_Spline_Gas[i], FD_Rat_Gas[i], nr_spline+1, 1e32, 1e32, D_Rat2_Gas[i]);
        SplineD(Rad_Spline_Gas[i], FSZ_Rat_Gas[i], nr_spline+1, 1e32, 1e32, SZ_Rat2_Gas[i]);
        
        //for (int j = 0; j < nr_spline; ++j)
        //{
        //    cout << "Spline   " << Rad_Spline[j] << " " << FD_Rat[j] << " " 
        //         << D_Rat2[i] << " " << FSZ_Rat[j] << " " << SZ_Rat2[i] << endl;
        //}
        
        // KK: Make niter iterations on the correction functions for 
        // midplane density and vertical velocity dispersion.
        // DP: Integration is compared to model disk structure from the density
        // functions. The correction functions measure the deviation from the
        // target disk density given by DiskDensfI().
        for (int k = 0; k < n_iter; ++k)
        {
            for (int j = 0; j < nr_spline+1; ++j)
            {
                double rad = Rad_Spline_Gas[i][j];
                
                GetOmegaKappa(rad, freq_omega, freq_kappa);
                //cout << "freqs   " << k << " " << j << " " << freq_omega << " " 
                //     << freq_kappa << endl;
                
                double sig_r = sqrt(SigR2Gas(rad,i));
                
                double v_circ = rad*freq_omega;
                
                double d0 = 0, dz2 = 0;
                double dvr = 0.1*sig_r;
        
                // KK: integrate analytically over vr and vz, numerically over vphi
                for (int ivt = 1; ivt < 102; ++ivt)
                {
                    double vt = v_circ+(ivt-51)*dvr;
                    double df = DiskDF5intezGas(vt, rad, 0, i);
                    double c = 0.33333333333333*((ivt%2)*2+2);
                    d0 += c*dvr*df;
                    double df1 = DiskDF5intezGas(vt, rad, G.Z_GasDisk[i], i);
                    dz2 += c*dvr*df1;
                }
                
                //if (d0==0 && dz2==0)//because this happens when exp(arg) falls below machine limitations
                //{
                //    d0 = 1e-280;
                //    dz2 = 0.419974e-280;
                //}
                
                double dens_const = GetDensConst(rad, i);
                
                D_Rat_Gas[i][j] = d0/GasDiskDensfI2(rad, 0, i);
                DZ2_Rat_Gas[i][j] = dz2/GasDiskDensfI2(rad, G.Z_GasDisk[i], i);
                FZ_Rat_Gas[i][j] = log(GasDiskDensfI2(rad, 0, i)/
                                       GasDiskDensfI2(rad, G.Z_GasDisk[i], i))/log(d0/dz2);
            }
            
            //Output a table of residuals to the console
            cout << endl;
            
            for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
            {
                cout << setprecision(5) << left << Rad_Spline_Gas[i][j] << " ";
            }
            
            cout << endl;
            cout << "0.000 ";
            
            for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
            {
                cout << setprecision(5) << left << D_Rat_Gas[i][j] << " ";
            }
            
            cout << endl;
            cout << G.Z_GasDisk[i] << " ";
            
            for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
            {
                cout << setprecision(5) << left << DZ2_Rat_Gas[i][j] << " ";
            }
            
            cout << endl;
            
            for (int j = 0; j < nr_spline+1; ++j)
            {
                FD_Rat_Gas[i][j] = FD_Rat_Gas[i][j]/D_Rat_Gas[i][j];
                FSZ_Rat_Gas[i][j] = FSZ_Rat_Gas[i][j]/FZ_Rat_Gas[i][j];
            }
            
            SplineD(Rad_Spline_Gas[i], FD_Rat_Gas[i], nr_spline+1, 1e32, 1e32, D_Rat2_Gas[i]);
            SplineD(Rad_Spline_Gas[i], FSZ_Rat_Gas[i], nr_spline+1, 1e32, 1e32, SZ_Rat2_Gas[i]);
        }
        
        // KK: Finally, write out a more detailed table of residuals.
        for (int j = 0; j < nr_spline+1; ++j)
        {
            double rad = Rad_Spline_Gas[i][j];

            GetOmegaKappa(rad, freq_omega, freq_kappa);

            double sig_r = sqrt(SigR2Gas(rad,i));

            double v_circ = rad*freq_omega;

            double d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0;
            double dvr = 0.1*sig_r;

            for (int ivt = 1; ivt < 102; ++ivt)
            {
                double vt = v_circ+(ivt-51)*dvr;
                double df = DiskDF5intezGas(vt, rad, 0, i);
                double c = 0.33333333333333*((ivt%2)*2+2);
                d0 += c*dvr*df;
                df = DiskDF5intezGas(vt, rad, 0.5*G.Z_GasDisk[i], i);
                d1 += c*dvr*df;
                df = DiskDF5intezGas(vt, rad, G.Z_GasDisk[i], i);
                d2 += c*dvr*df;
                df = DiskDF5intezGas(vt, rad, 1.5*G.Z_GasDisk[i], i);
                d3 += c*dvr*df;
                df = DiskDF5intezGas(vt, rad, 2*G.Z_GasDisk[i], i);
                d4 += c*dvr*df;
            }
            
            double dens_const = GetDensConst(rad, i);
                
            dirat[0][j] = d0/GasDiskDensfI2(rad, 0, i);
            dirat[1][j] = d1/GasDiskDensfI2(rad, 0.5*G.Z_GasDisk[i], i);
            dirat[2][j] = d2/GasDiskDensfI2(rad, G.Z_GasDisk[i], i);
            dirat[3][j] = d3/GasDiskDensfI2(rad, 1.5*G.Z_GasDisk[i], i);
            dirat[4][j] = d4/GasDiskDensfI2(rad, 2*G.Z_GasDisk[i], i);
        }
        
        //Output a table of residuals to the console
        cout << endl;
        
        for (int j = 0; j < nr_spline; j+=nr_spline/10)
        {
            cout << setprecision(5) << left << Rad_Spline_Gas[i][j] << " ";   
        }
        
        cout << endl;
        for (int k = 0; k < 5; ++k)
        {
            cout << 0.5*k*G.Z_GasDisk[i] << " ";
        
            for (int j = 0; j < nr_spline+1; j+=nr_spline/10)
            {
                cout << setprecision(5) << left << dirat[k][j] << " ";
            }
        
            cout << endl;
        }
        
        cout << endl;
            
        //write to cordbh.dat for each disk component if requested
        if (do_file_io)
        {
            cordbhout << G.Sigma_0_Gas[i] << " " << G.R_Sigma_Gas[i] << " " << nr_spline << endl;
            
            for (int j = 0; j < nr_spline+1; ++j)
            {
                cordbhout << Rad_Spline_Gas[i][j] << " " << D_Rat_Gas[i][j] << " " 
                          << D_Rat2_Gas[i][j] << " " << FD_Rat_Gas[i][j] << " " 
                          << FZ_Rat_Gas[i][j] << " " << FSZ_Rat_Gas[i][j] << endl;
            }
        }
        else
        {
            //for (int j = 0; j < nr_spline+1; ++j)
            //{
            //    cout << setw(12) << Rad_Spline[i][j] << " " << setw(12) 
            //         << D_Rat[i][j] << " " << setw(12) << D_Rat2[i][j] << " "
            //         << setw(12) << FD_Rat[i][j] << " " << setw(12) 
            //         << FZ_Rat[i][j] << " " << setw(12) << FSZ_Rat[i][j] << endl;
            //}
        }
        
        //Now get the Toomre stability parameters Q and X
        double q_avg = 0, x_avg = 0, q_scale, x_scale;
        int points = 0;
        ofstream toomreq("toomre.dat"), toomrex("toomres_other_parameter.dat");
        
        for (int j = 0; j < nr_spline; ++j)
        {
            double rad = Rad_Spline_Gas[i][j];
            GetOmegaKappa(rad, freq_omega, freq_kappa);
            double sig_r = sqrt(SigR2Gas(rad,i));
            double sig_den = GasDiskDensfI2(rad,0,i)*2*G.Z_GasDisk[i];
            
            double mode = 2;
            double q = sig_r*freq_kappa/PI/sig_den;
            double x = freq_kappa*freq_kappa*rad*oneover2pi/mode/sig_den;
            
            if (rad < 2*G.R_GasDisk[i])
            {
                q_scale = q;
                x_scale = x;
            }
            
            if (rad > 0.5*G.R_GasDisk[i] && rad < min(3*G.R_GasDisk[i], GasDisk_Edge[i]))
            {
                q_avg += q;
                x_avg += x;
                ++points;
            }
            
            toomreq << rad << " " << q << endl;
            toomrex << rad << " " << x << endl;
        }
            
        Q_scale_gas.push_back(q_scale);
        X_scale_gas.push_back(x_scale);
        Q_avg_gas.push_back(q_avg/points);
        X_avg_gas.push_back(x_avg/points);
            
        cout << "Q at 2 scale lengths = " << Q_scale_gas[i] << endl;
        cout << "X at 2 scale lengths = " << X_scale_gas[i] << endl;
        cout << "Q average = " << Q_avg_gas[i] << endl;
        cout << "X average = " << X_avg_gas[i] << endl;
        cout << endl;
        
    }
    
    //Must also get a total Q and X
}
                                
        
    
    
