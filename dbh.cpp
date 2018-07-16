//The main function that yields the harmonics. To be called after all 
//parameters have been read with GetParameters() and memory allocated
//with AllocateVectors().

#include "galaxy.h"

void DBH(void)
{
    cout.precision(10);
    
    //Get Legendre constants
    for (int i = 0; i < 100; ++i)
        Plcon[i] = sqrt((2*i+1)*oneover4pi);
    
    //Get rho_0 for the halo
    if (halo_flag == 1)
    {
        halo_const = GetHaloConst();
    }
    
    //Get p, b, and rho_0 for the bulge
    if (bulge_flag == 1)
    {
        if (G.ppp < 0)
        {
            G.ppp = 1 - 0.6097/G.n_sersic + 0.05563/G.n_sersic/G.n_sersic;
            G.b_n = Get_b_n();
	        G.rho_0 = Get_rho_0();
            gamma_comp = gsl_sf_gamma(G.n_sersic*(3 - G.ppp));
        }
    }
    
    //Set the length scales by finding the max and min of all length scales 
    //and cutoff radii in the galaxy. dr is multipled by the smallest length
    //scale in the model and therefore becomes the smallest radius on the grid
    //while r_max becomes the largest radius.
    double min_disk_length = 1e20, min_halo_length = 1e20, min_bulge_length = 1e20;
    double max_disk_length = -1e20, max_halo_length = -1e20, max_bulge_length = -1e20;
    
    if (disk_flag == 1)
    {
        for (int i = 0; i < G.R_Disk.size(); ++i)
        {
            if (G.R_Disk[i] < min_disk_length)
                min_disk_length = G.R_Disk[i];
            if (G.Out_Disk[i]+5*G.Dr_Trunc[i] > max_disk_length)
                max_disk_length = G.Out_Disk[i]+5*G.Dr_Trunc[i];            
        }
    }
    
    if (halo_flag == 1)
    {
        if (G.a_halo < min_halo_length)
            min_halo_length = G.a_halo;
        if (G.c_halo+5*G.drtrunc_halo > max_halo_length)
            max_halo_length = G.c_halo+5*G.drtrunc_halo;
    }
    
    if (bulge_flag == 1)
    {
        if (G.a_bulge < min_bulge_length)
            min_bulge_length = G.a_bulge;
        if (G.a_bulge > max_bulge_length)
            max_bulge_length = G.a_bulge;
    }   
    
    //dr *= min(min_halo_length, min(min_bulge_length, min_disk_length));
    r_max = 2*max(max_halo_length, max(max_bulge_length, max_disk_length));
    log_rmax = log10(r_max);//to account for vector overflows right at the edge
    log_dr = log10(dr);
    delta_logr = (log_rmax - log_dr) / (nr-2);
    dr=r_max/(nr-1);
    
#ifdef DEBUG    
    cout << "Set length scale " << dr << " " << r_max << " " << nr << " " 
         << max_halo_length << " " << max_disk_length << " " << max_bulge_length << endl;
    cout << "Set length scale " << dr << " " << r_max << " " << nr << " " 
         << min_halo_length << " " << min_disk_length << " " << min_bulge_length << endl;
    cout << "Input parameters: " << endl;
    if (halo_flag)
    {
        cout << "  Halo: " << G.c_halo << " " << G.v_halo << " " << G.a_halo << " "
             << G.drtrunc_halo << " " << G.cusp << " " << G.halo_stream << endl;
    }
    if (disk_flag)
    {
        cout << "  Disk: " << G.M_Disk[0] << " " << G.R_Disk[0] << " " << G.Out_Disk[0]
             << " " << G.Z_Disk[0] << " " << G.Dr_Trunc[0] << " " << G.Sigma_0[0] 
             << " " << G.R_Sigma[0] << endl;
    }
    if (bulge_flag)
    {
        cout << "  Bulge: " << G.n_sersic << " " << G.ppp << " " << G.v_bulge << " "
             << G.a_bulge << " " << G.bulge_stream << endl; 
    }
#endif
    
    //Set radius vector
    Radius[0] = 0.0;
    
    for(int i = 1; i < nr; ++i)
    {
        double rad = (i-1)*delta_logr + log_dr;
        Radius[i] = i*dr;//pow(10,rad);
        //cout<<"radius "<<i<<" "<<Radius[i]<<" "<<r_max<<endl;
    }
    
#ifdef DEBUG    
    cout << "Set radius and length scale" << endl;
#endif
    
    CopyGlobalsToDevice1();
    
    //Make a bunch of function calls for the Eddington inversion
    if (halo_flag) 
    {
        if (contraction_flag)
        {
            ContractHalo();
            HaloMassDwithRadius();
        }
        
        cout << "    Calculated initial estimate of halo potential" << endl;        
        HaloPotentialEstimate();
        cout << "    Calculated initial estimate of halo potential" << endl;        
    }
    
    //CopyGlobalsToDevice1();
    if (disk_flag) 
    {
        DiskPotentialEstimate();
        cout << "    Calculated initial estimate of disk potential" << endl;        
    }
    
    if (gasdisk_flag) 
    {
        GasDiskPotentialEstimate();
        cout << "    Calculated initial estimate of gas disk potential" << endl;        
    }
    
    if (bulge_flag) 
    {
        BulgePotentialEstimate();
        cout << "    Calculated initial estimate of bulge potential" << endl;        
    }
    
    GenTableE();
    cout << "    Constructed table of energies" << endl;
    CopyGlobalsToDevice2();
    
    if (bulge_flag) 
    {
        //GenSersicDistFunc();
        GenSersicDistFuncCUDA();
        cout << "    Calculated Sersic distribution function" << endl;
        //exit(0);       
    }
    
    if (halo_flag == 1) 
    {
        GenNFWDistFunc();
        //GenNFWDistFuncCUDA();
        cout << "    Calculated halo distribution function" << endl;
    }

    if (halo_flag == 1) 
    {
        GenDensPsiHalo();    
        cout << "    Calculated halo density from energy" << endl;
    }
    
    if (bulge_flag == 1) 
    {
        GenDensPsiBulge();
        cout << "    Calculated bulge density from energy" << endl;
    }    
    
    for(int ir = 0; ir < nr; ++ir)
    {
        double rad = Radius[ir];
        A_Pot[0][ir] = sqrt4pi * GetTotalPsi(rad);
        //cout<<ir<<" "<<A_Pot[0][ir] << endl;
        
        for (int l = 2; l < l_max+1; l+=2)
            A_Pot[l/2][ir] = 0;
    }
    
    if (gasdisk_flag)
    {
        GasDiskCUDA();
        //CalculatePolytropeConstants();
        //CalculateGasDensConstants();
        //CalculateGasScaleHeights();
        //GetAppGasDiskPot();
        for (int i = 0; i < nr; i+=1000)
        {
            cout << "Did " << i << " calculation " << setw(12) << Radius[i]
                 << " " << GasDensity_Const[0][i] << "               "
                 << Pot(Radius[i], 0) << " " << endl;
        }
    }
    
    /***********************************
    double masst = 0;
    for (int i = 1; i < nr_ac; ++i)
    {
        double rad = Halo_Original_Radius[i], z = 0;
        double radrad = sqrt(rad*rad+z*z);
        double rrat = radrad/G.a_bulge;
        double targetdens = G.rho_0*pow(rrat, -G.ppp)*exp(-G.b_n*pow(rrat, 1.0/G.n_sersic));
        
        double volume = fourpi/3*(pow(rad, 3)-pow(Halo_Original_Radius[i-1], 3));
        masst += BulgeDens(rad, z)*volume;
        
        cout << rad/r_max << "  " << targetdens << " "
             << BulgeDens(rad, z) << " " << masst/2466.246138 << "    " 
             << GetBulgeMass(radrad)/2466.246138 << endl;
    }
    ***********************************/
    
    //Note that l_max and lmax are different. l_max is the overall maximum l
    //while lmax is the nested maximum l that gets iterated up to l_max
    
    int lmaxstep = 2, lmax = 0;
    int n_theta = max(10, lmax*10+2);
    
    //KK: now iterate. number of iterations and change in lmax depends on
    //initial conditions.  iterate on first cycle until tidal radius is
    //stable, then once for every harmonic added, and until convergence
    //once l has reached the desired maximum lmax.
      
    double drtidal = 2*dr, rtidal_old = 1e30, rtidal, r_old;
    double tidal_check = rtidal_old;
    int lmax_old = -2, iter_outside = 0;

    cout << "    Starting main loop" << endl;
    
/******************************************************************************/
/******************************MAIN LOOP***************************************/
/******************************************************************************/
    
    for (iter = 0; iter < 1000; ++iter)
    {        
        if (gasdisk_flag && iter > 0)
        {
            GasDiskCUDA();
            //CalculatePolytropeConstants();
            //CalculateGasDensConstants();
            //CalculateGasScaleHeights();
            //GetAppGasDiskPot();
            cout << "    Calculated gas constants" << endl;
            for (int i = 0; i < nr; i+=1000)
            {
                cout << "Did " << i << " calculation " << setw(12) << Radius[i]
                     << " " << GasDensity_Const[0][i] << "               "
                     << Pot(Radius[i], 0) << " " << endl;
            }
        }
        //GetAppDiskPot();
        
        ////////////////////////////////////////////////////
        //for (int i = 0; i < nr; i+=100)
        //{
        //    cout << Radius[i] << " " << Pot(0, Radius[i]) << endl;
        //}
           
        if (lmax == 0 || lmax == l_max)
        {
            if (drtidal < dr && iter > 10)
            {
                lmax = lmax + lmaxstep;
                n_theta = max(10, lmax*5);//lmax*4 + 2;
            }
        }
        else
        {
            lmax = lmax + lmaxstep;
            n_theta = max(10, lmax*5);//lmax*4 + 2;
        }
        
        if (lmax == l_max+lmaxstep)
            break;
        
        //KK: Now get the harmonics of the density in this potential --> A_Dens
        //NB that dens only gives the density without an approximate sech**2
        //component --- that part of the potential is not represented by the
        //harmonics.  The function dens(r,z) returns the density -
        //high-frequency cpt The full density is returned by totdens
        
        double eps = 0.00001;
        A_Dens[0][0] = Dens(eps, 0) * sqrt4pi;
        //cout<<"adens   "<<A_Dens[0][0]<<endl;
        
        for (int l = 2; l < lmax+1; l+=2)
            A_Dens[l/2][0] = 0;
        
        for (int ir = 1; ir < nr; ++ir)
            A_Dens[0][ir] = 0;
        
        //KK: nr_limit will mark the outermost radial bin with non-zero density.
        nr_limit = nr;
        
        //KK: integrate density * spherical harmonic function over quadrant use
        //cos(theta) as independent variable.
        //DP: d_cos_theta is d(cos(theta))
        
        for (int l = 0; l < lmax+1; l+=2)
        {
            GetADensCUDA(l);
                        
            //for (int ir = 1; ir < nr_limit; ++ir)
            //{
                
                //KK: mark the first even radial bin on which the density has fallen to zero.
                //if (l == 0 && s == 0)
                //{
                //    nr_limit = nr;
                //    nr_limit += nr_limit%2;
                //    break;
                //}
            //}            
        }
        
        //cout<<"apot00   "<<A_Pot[0][0]<<" "<<A_Pot[0][2]<<" "<<A_Pot[0][4] 
        //    << "   "<<A_Dens[0][0]<<" "<<A_Dens[0][2]<<" "<<A_Dens[0][4]<<endl;    
        
        //KK: now get the potential harmonics of this new density. (BT 2-208)
        //Simpson's rule integration.
        
        //cout << "Done adens" << endl;
                
        for (int l = 0; l < lmax+1; l+=2)
        {
            S_1[0] = 0;
            double rad = Radius[2];
            S_1[2] = (rad*dr/3)*(4*A_Dens[l/2][1]*pow(1-dr/rad,l+2) + A_Dens[l/2][2]);
            r_old = rad;
            
            for (int ir = 4; ir < nr; ir+=2)
            {
                double rad = Radius[ir];  
                double s1a = (rad*dr/3)*(A_Dens[l/2][ir-2]*pow(1-2*dr/rad,l+2) +
                             4*A_Dens[l/2][ir-1]*pow(1-dr/rad,l+2) + A_Dens[l/2][ir]);
                S_1[ir] = s1a + S_1[ir-2]*pow(r_old/rad, l+1);
                r_old = rad;
            }
            
            S_2[nr-1] = 0;
            r_old = Radius[nr-1];
            
            for (int ir = nr-3; ir > 1; ir-=2)
            {
                double rad = Radius[ir];  
                double s2a = (rad*dr/3) * (A_Dens[l/2][ir+2]*pow(1+2*dr/rad,1-l) +
                             4*A_Dens[l/2][ir+1]*pow(1+dr/rad,1-l) + A_Dens[l/2][ir]);
                S_2[ir] = s2a + S_2[ir+2]*pow(rad/r_old, l);
                r_old = rad;
            }
         
            //KK: replace the potential harmonics with a mean of the previous
            //iteration (25%) and the current one (75%). This damps out
            //oscillations that otherwise occur.  if this is the first time this
            //harmonic is calculated, use the entire new value.
            
            //DP: The 25% and 75% figures above are in fact reversed in the
            //expression below. But it still works.

            for (int ir = 2; ir < nr; ir+=2)
            {
                if (l <= lmax_old)
                {
                    A_Pot[l/2][ir] = fraction*A_Pot[l/2][ir] + 
                                     (1-fraction)*fourpi*(S_1[ir]+S_2[ir])/(2.0*l+1);
                }
                else
                {
                    A_Pot[l/2][ir] = fourpi*(S_1[ir]+S_2[ir])/(2.0*l+1);
                }
            }
            
            //KK: Calculate the 1st and 2nd-order radial gradients
            
            for (int ir = 2; ir < nr; ir+=2)
            {
                double rad = Radius[ir];
                F_R[l/2][ir] = -fourpi*(-(l+1)*S_1[ir]+l*S_2[ir]) / (2.0*l+1) / rad;
                FR2[l/2][ir] = -fourpi*((l+1)*(l+2)*S_1[ir]/rad/rad + l*(l-1)*S_2[ir]/rad/rad - 
                               (2.0*l+1)*A_Dens[l/2][ir])/(2.0*l+1);
            }
        }
        
        // KK: now interpolate the gaps first quadratically interpolate the
        // monopole back to the origin.  the remaining multipoles are zero
        // there.
        
        // DP: GSL polynomial spline algorithm
        
//         int n_poly = 5;
//         double r_zero = 0;
//         double xa[5] = {2, 4, 6, 8, 10};
//         double ya[5] = {A_Pot[0][2], A_Pot[0][4], A_Pot[0][6], A_Pot[0][8], A_Pot[0][10]};
//         double dy;
//         
//         gsl_interp_accel *acc = gsl_interp_accel_alloc();
//         const gsl_interp_type *T = gsl_interp_polynomial;
//         gsl_spline *spline = gsl_spline_alloc(T, n_poly);
//         gsl_spline_init(spline, xa, ya, n_poly);
//         
//         A_Pot[1][0] = gsl_spline_eval(spline, r_zero, acc);
//         
//         gsl_spline_free(spline);
//         gsl_interp_accel_free(acc);
        
        A_Pot[0][0] = 3*(A_Pot[0][2]-A_Pot[0][4])+A_Pot[0][6];
        
        FR2[0][0] = 2*FR2[0][2] - FR2[0][4];
        
        for (int l = 2; l < lmax+1; l+=2)
        {
            A_Pot[l/2][0] = 0;
            F_R[l/2][0] = 0;
            FR2[l/2][0] = 0;
        }
        
        //KK: then linearly interpolate other bins.
        
        for (int ir = 1; ir < nr; ir+=2)
        {
            for (int l = 0; l < lmax+1; l+=2)
            {
                A_Pot[l/2][ir] = (A_Pot[l/2][ir-1] + A_Pot[l/2][ir+1]) / 2;
                F_R[l/2][ir] = (F_R[l/2][ir-1] + F_R[l/2][ir+1]) / 2;             
                FR2[l/2][ir] = (FR2[l/2][ir-1] + FR2[l/2][ir+1]) / 2;
                //if(l>0)
                //cout << "harmonics "<<l<<" "<<dr*ir<<" "<<A_Pot[l/2][ir]
                //     <<" "<<F_R[l/2][ir]<<" "<<FR2[l/2][ir]<<endl;
            }
        }
                
        //KK: finally reset the potential at the origin to psi0
        //Note that the fake disk potential is zero at the origin.
        if (halo_flag || bulge_flag)
        {
            double a00 = A_Pot[0][0];
            //cout<<a00<<" "<<psi_crit<<" "<<psi_0<<endl;
            
            for (int ir = 0; ir < nr; ++ir)
            {
                A_Pot[0][ir] = A_Pot[0][ir] + psi_0*sqrt4pi - a00;
                /* if (iter>0) */ 
                //     cout<<"apot  "<<ir<<" "<<A_Pot[0][ir]<<" "<<psi_0*sqrt4pi <<" "<< a00<<endl;
            }
            
            if (a00*oneoversqrt4pi - psi_0 > psi_crit)
            {
                cout << "Iter " << iter << ":   lmax = " << lmax 
                     << "   tidal radius = infinite " << psi_crit << " "
                     << psi_0 << " " << a00 << endl;
                ++iter_outside;
                
                drtidal = 2*dr;
            }
            else
            {
                double pot_r = psi_0;
                
                //KK: look for new tidal radius of the model defined as the radius at
                //the equator where the potential is equal to psic
      
                for (int ir = 1; ir < nr; ++ir)
                {
                    double potrm1 = pot_r;
                    double rad = Radius[ir];
                    pot_r = Pot(rad, 0);
                
                    double aa = potrm1 - psi_crit;
                    double bb = pot_r - psi_crit;
                
                    if (aa*bb <= 0)
                    {
                        double dpot = pot_r - potrm1;
                        
                        if (dpot == 0)
                        {
                            rtidal = Radius[ir-1];
                        }
                        else
                        {
                            rtidal = (ir - 1 - aa/dpot)*dr;
                        }
                    
                        drtidal = fabs(rtidal - rtidal_old);
                        tidal_check = fabs(rtidal - rtidal_old)/rtidal;
                    
                        cout << "Iter " << iter << ":   lmax = " << lmax 
                             << "   tidal radius = " << rtidal << " " << psi_crit << " "
                             << psi_0 << " " << a00 << endl;
                    
//                         for (int k = 0; k < nr-1; k+=100)
//                         {
//                             double rad = Radius[k];
//                             cout << "     " << dr << " " << rad << " " << Pot(rad, 0) << " "
//                                  << Pot(rad, dr) << " " << Pot(rad, 10*dr) << " " 
//                                  << Pot(rad, 100*dr) << endl;
//                         } 
                        
                        rtidal_old = rtidal;
                    
                        break;
                    }
                    else if (ir == nr-1)
                    {
                        cout << "Iter " << iter << ":   lmax = " << lmax 
                             << "   tidal radius = outside grid " << psi_crit << " "
                             << psi_0 << " " << a00 << endl;    
                               
                        drtidal = 2*dr;
                        break;
                    }
                }
            }
        }
        
                
        //KK: write out the changes in the potential at the reference points 
        //at this iteration.
        //DP: Apparently refpoints is never used so this is commented out
        //if (disk_flag == 1)
        //{
        //    for (int iref = 0; iref < 10; ++iref)
        //    {
        //        old_pref[iref] = pref[iref];
        //        pref[iref] = pot(rref[iref], zref[iref]);
        //        
        //        refpoints << iter << " " << lmax << " " 
        //                  << pref(iref) - oldpref(iref) << endl;
        //    }
        //}
           
        //KK: now repeat with this new potential!
        lmax_old = lmax;
        
        //cout << "Finished iteration" << endl;
    }    
                
/*********************************************************************************/
/*****************************END OF MAIN LOOP************************************/
/*********************************************************************************/

    cout << "diskdens " << DiskDensf(0,0) << endl;iter++;
        
    cout << nr << endl;                    
    total_mass = F_R[0][nr-1]*Radius[nr-1]*Radius[nr-1]*oneoversqrt4pi;
    cout << "Total mass = " << total_mass << endl;
    
    //KK: Calculate force and potential for halo only
    
    halo_mass = 0, halo_edge = 0, bulge_mass = 0, bulge_edge = 0;
    disk_mass = 0, disk_edge = 0;
    lmax = l_max;
        
    if (halo_flag)
        HaloPotential();
    if (bulge_flag)
        BulgePotential();
    if (smbh_flag)
        GenBlackHole(G.bh_mass);
    
    if (disk_flag)
    {
        total_disk_mass = total_mass - halo_mass - bulge_mass;
        
        for (int i = 0; i < disk; ++i)
        {
            Disk_Edge.push_back(G.Out_Disk[i] + 2*G.Dr_Trunc[i]);
            Disk_Mass.push_back(GetDiskMass(i));
            
            cout << "Disk " << i+1 << " mass = " << Disk_Mass[i] << endl;
            cout << "Disk " << i+1 << " edge = " << Disk_Edge[i] << endl;
            
            if (disk_edge < Disk_Edge[i])
            {
                disk_edge = Disk_Edge[i];
            }
        }
        
        cout << "Total disk mass = " << total_disk_mass << endl;
        cout << "Outermost disk radius = " << disk_edge << endl;
    }
      
    if (gasdisk_flag)
    {
        total_disk_mass = total_mass - halo_mass - bulge_mass;
        
        for (int i = 0; i < gas_disk; ++i)
        {
            GasDisk_Edge.push_back(G.Out_GasDisk[i] + 2*G.Dr_Trunc_Gas[i]);
            GasDisk_Mass.push_back(GetGasDiskMass(i));
            
            cout << "GasDisk " << i+1 << " mass = " << GasDisk_Mass[i] << endl;
            cout << "GasDisk " << i+1 << " edge = " << GasDisk_Edge[i] << endl;
            
            if (disk_edge < GasDisk_Edge[i])
            {
                disk_edge = GasDisk_Edge[i];
            }
        }
        
        cout << "Total disk mass = " << total_disk_mass << endl;
        cout << "Outermost disk radius = " << disk_edge << endl;
    }
      
    //DP: Output to mr.dat the masses and radii and to dbh.dat the harmonics
    //if requested
    
    if (do_file_io)
    {
	    ofstream tidal_file("rtidal.dat", ios::out);
		tidal_file << rtidal << endl;
        WriteDBHDat();
    }
    
//     double mass = 0, massh = 0;
//     for (int i = 0; i < nr; ++i)//i+=10)
//     {
//         double r = Radius[i];
//         //double poly = GetPolyConst(r,0);
//         double poly = 0;//GetDensConst(r,0);
//         double surfden = 0;//GasDiskSurfaceDensfI2(r, 0, poly);
//         //double haloden = HaloDens(r,0);
//         //mass += r*dr*2*PI*surfden;
//         
//         double masshh = 0;
//         for (int j = 1; j < nr; ++j)
//         {
//             double zz = j*dr;
//             masshh += dr*HaloDens(r, zz);
//             surfden += dr*GasDiskDensf2(r, zz);
//         }
//         massh += 2*PI*dr*r*masshh;
//         massh += PI*dr*r*dr*HaloDens(r, 0);
//         mass += 2*PI*dr*r*surfden;
//         mass += PI*dr*r*dr*GasDiskDensf2(r, 0);
//         surfden += 0.5*dr*GasDiskDensf2(r, 0);
// 
//         cout << Radius[i] << " " << 2*surfden 
//              << " " <<  GasDisk_Const[0]*exp(-r/G.R_GasDisk[0]) << "    "
//              << 2*mass << " " << 2*massh << endl;
//         
// //         for (double z = 0; z < 10; z+=dr)
// //         {
// //             double psi = Pot(r,z);
// //             double partialsurfden = GasDiskSurfaceDensfI2(r, 0, poly, z);
// //             double ratio = partialsurfden/surfden;
// //             
// //             //if(ratio>0.68)
// //             {
// //                 cout << "    " << r << " " << setw(12) << z << " "
// //                      << GasDiskDensINoTrunc2(r,z,psi,0,poly) << " "
// //                      << partialsurfden << " " << ratio << "       "
// //                      << DiskDens(r,z,psi) << endl;
// //                 //break;
// //             }
// //         }
//     }
//     
//     exit(0);
    
    //Check to see if Poisson's equation is satisfied
    //First tabulate the partial derivatives
//     int fac = 30, starter = 0;
//     vector<vector<vector<double> > > PsiDerivX, PsiDerivY, PsiDerivZ;//[nr/fac][nr/fac][nr/fac];
//     PsiDerivX.resize(nr/fac);
//     PsiDerivY.resize(nr/fac);
//     PsiDerivZ.resize(nr/fac);
//     
//     for(int i=0; i < PsiDerivX.size(); ++i)
//     {
//         PsiDerivX[i].resize(nr/fac, vector<double>(nr/fac, 0));
//         PsiDerivY[i].resize(nr/fac, vector<double>(nr/fac, 0));
//         PsiDerivZ[i].resize(nr/fac, vector<double>(nr/fac, 0));
//     }
//     
//     cout << "Got vectors" << endl;
//     
//     for (int i = 1; i < nr/fac; ++i)
//     {
//         for (int j = 1; j < nr/fac; ++j)
//         {
//             for (int k = 1; k < nr/fac; ++k)
//             {
//                 double x = dr*(i+starter), y = dr*(j), z = dr*(k);
//                 double r = sqrt(x*x+y*y);
//                 //double rhs = fourpi*TotalDens(r, z);
//                 
//                 double delta = dr;
//                 double x1 = x+delta, x2 = x-delta, dx = 2*delta; 
//                 double y1 = y+delta, y2 = y-delta, dy = 2*delta; 
//                 double z1 = z+delta, z2 = z-delta, dz = 2*delta; 
//                 
//                 double rx1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 double ry1 = sqrt(x*x+y1*y1), ry2 = sqrt(x*x+y2*y2);
//                 //double rz1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 
//                 //Get dpsi/dx
//                 double psix1 = Pot(rx1, z), psix2 = Pot(rx2, z);
//                 double psiy1 = Pot(ry1, z), psiy2 = Pot(ry2, z);
//                 double psiz1 = Pot(r, z1), psiz2 = Pot(r, z2);
//                 
//                 double dpsi_dx = (psix1-psix2)/(x1-x2);
//                 double dpsi_dy = (psiy1-psiy2)/(y1-y2);
//                 double dpsi_dz = (psiz1-psiz2)/(z1-z2);
//                 
//                 PsiDerivX[i][j][k] = dpsi_dx;
//                 PsiDerivY[i][j][k] = dpsi_dy;
//                 PsiDerivZ[i][j][k] = dpsi_dz;
//                 
//                 if (k%100==0&&j%100==0&&i%10==0) cout << i << " " << j << " " << " " << k << endl;
//                 //cout << i << " " << j << " " << " " << k << "    "
//                 //     << dpsi_dx << " " << dpsi_dy << " " << dpsi_dz
//                 //     << endl;
//             }
//         }
//     }
//     
//     cout << "Part 1 done" << endl;
//     
//     //Now get the Laplacian and compare to the RHS of Poisson's equation
//     for (int i = 2; i < nr/fac-1; ++i)
//     {
//         for (int k = 2; k < nr/fac-1; ++k)
//         {
//             for (int j = 2; j < nr/fac-1; ++j)
//             {
//                 double x = dr*(i+starter), y = dr*(j), z = dr*(k);
//                 double r = sqrt(x*x+y*y);
//                 double rhs = fourpi*TotalDens(r, z);
//                 
//                 double delta = dr;
//                 double x1 = x+delta, x2 = x-delta, dx = 2*delta; 
//                 double y1 = y+delta, y2 = y-delta, dy = 2*delta; 
//                 double z1 = z+delta, z2 = z-delta, dz = 2*delta; 
//                 
//                 //double rx1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 //double ry1 = sqrt(x*x+y1*y1), ry2 = sqrt(x*x+y2*y2);
//                 //double rz1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 
//                 //Get dpsi/dx
//                 double dpsidx1 = PsiDerivX[i+1][j][k], dpsidx2 = PsiDerivX[i-1][j][k];
//                 double dpsidy1 = PsiDerivY[i][j+1][k], dpsidy2 = PsiDerivY[i][j-1][k];
//                 double dpsidz1 = PsiDerivZ[i][j][k+1], dpsidz2 = PsiDerivZ[i][j][k-1];
//                 
//                 double d2psi_dx2 = (dpsidx1-dpsidx2)/(x1-x2);
//                 double d2psi_dy2 = (dpsidy1-dpsidy2)/(y1-y2);
//                 double d2psi_dz2 = (dpsidz1-dpsidz2)/(z1-z2);
//                 
//                 double laplacian = d2psi_dx2 + d2psi_dy2 + d2psi_dz2;
//                 double rrad = sqrt(r*r+z*z);
//                 
//                 //if (z>0 && z<0.2 && x<0.5 && x>-0.3)
//                 if ((z>1 && z<1.1 && x<0.4 && x>0.3) || 
//                     (z>0 && z<0.1 && x<0.4 && x>0.3))
//                 cout << x << " " << y << " " << " " << z << "    "
//                      << d2psi_dx2 << " " << d2psi_dy2 << " " << d2psi_dz2
//                      << "     " << -laplacian << " " << rhs << "       " 
//                      << fourpi*DiskDensf(r, z) << " " 
//                      << fourpi*HaloDens(r, z) << " "
//                      << Pot(r, z) << endl;
//                      //<< fourpi*RawHaloProfile(rrad) << endl;
//             }
//         }
//     }
//     
//     cout << "Part 2 done" << endl;
//     
//     //Now get the Laplacian and compare to the RHS of Poisson's equation
//     for (int i = 2; i < nr/fac-1; ++i)
//     {
//         for (int j = 2; j < nr/fac-1; ++j)
//         {
//             for (int k = 2; k < nr/fac-1; ++k)
//             {
//                 double x = dr*(i+starter), y = dr*(j), z = dr*(k);
//                 double r = sqrt(x*x+y*y);
//                 double rhs = fourpi*TotalDens(r, z);
//                 
//                 double delta = dr;
//                 double x1 = x+delta, x2 = x-delta, dx = 2*delta; 
//                 double y1 = y+delta, y2 = y-delta, dy = 2*delta; 
//                 double z1 = z+delta, z2 = z-delta, dz = 2*delta; 
//                 
//                 //double rx1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 //double ry1 = sqrt(x*x+y1*y1), ry2 = sqrt(x*x+y2*y2);
//                 //double rz1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 
//                 //Get dpsi/dx
//                 double dpsidx1 = PsiDerivX[i+1][j][k], dpsidx2 = PsiDerivX[i-1][j][k];
//                 double dpsidy1 = PsiDerivY[i][j+1][k], dpsidy2 = PsiDerivY[i][j-1][k];
//                 double dpsidz1 = PsiDerivZ[i][j][k+1], dpsidz2 = PsiDerivZ[i][j][k-1];
//                 
//                 double d2psi_dx2 = (dpsidx1-dpsidx2)/(x1-x2);
//                 double d2psi_dy2 = (dpsidy1-dpsidy2)/(y1-y2);
//                 double d2psi_dz2 = (dpsidz1-dpsidz2)/(z1-z2);
//                 
//                 double laplacian = d2psi_dx2 + d2psi_dy2 + d2psi_dz2;
//                 double rrad = sqrt(r*r+z*z);
//                 
//                 //if (y>0 && y<2 && x<0.5 && x>-0.3)
//                 if (y>0 && y<0.1 && x<0.5 && x>0.3)
//                 cout << x << " " << y << " " << " " << z << "    "
//                      << d2psi_dx2 << " " << d2psi_dy2 << " " << d2psi_dz2
//                      << "     " << -laplacian << " " << rhs << "       " 
//                      << fourpi*DiskDensf(r, z) << " "
//                      << fourpi*HaloDens(r, z) << " "
//                      << Pot(r, z)  << endl;
//                      //<< fourpi*RawHaloProfile(rrad) << endl;
//             }
//         }
//     }
    
//     starter = 250;
//     for (int i = 2; i < nr/fac-1; ++i)
//     {
//         for (int j = 2; j < nr/fac-1; ++j)
//         {
//             for (int k = 2; k < nr/fac-1; ++k)
//             {
//                 double x = dr*(i), y = dr*(j), z = dr*(k+starter);
//                 double r = sqrt(x*x+y*y);
//                 double rhs = fourpi*TotalDens(r, z);
//                 
//                 double delta = dr;
//                 double x1 = x+delta, x2 = x-delta, dx = 2*delta; 
//                 double y1 = y+delta, y2 = y-delta, dy = 2*delta; 
//                 double z1 = z+delta, z2 = z-delta, dz = 2*delta; 
//                 
//                 //double rx1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 //double ry1 = sqrt(x*x+y1*y1), ry2 = sqrt(x*x+y2*y2);
//                 //double rz1 = sqrt(x1*x1+y*y), rx2 = sqrt(x2*x2+y*y);
//                 
//                 //Get dpsi/dx
//                 double dpsidx1 = PsiDerivX[i+1][j][k], dpsidx2 = PsiDerivX[i-1][j][k];
//                 double dpsidy1 = PsiDerivY[i][j+1][k], dpsidy2 = PsiDerivY[i][j-1][k];
//                 double dpsidz1 = PsiDerivZ[i][j][k+1], dpsidz2 = PsiDerivZ[i][j][k-1];
//                 
//                 double d2psi_dx2 = (dpsidx1-dpsidx2)/(x1-x2);
//                 double d2psi_dy2 = (dpsidy1-dpsidy2)/(y1-y2);
//                 double d2psi_dz2 = (dpsidz1-dpsidz2)/(z1-z2);
//                 
//                 double laplacian = d2psi_dx2 + d2psi_dy2 + d2psi_dz2;
//                 
//                 if (x < 0.8 && x > 0.7)
//                 cout << x << " " << y << " " << " " << z << "    "
//                      << d2psi_dx2 << " " << d2psi_dy2 << " " << d2psi_dz2
//                      << "     " << -laplacian << " " << rhs << endl;
//             }
//         }
//     }
}
     
     
     
