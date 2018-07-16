#include "galaxy.h"

void ReadDensPsiBulge(void)
{
    ifstream denspsibulge("denspsibulge.dat");
    ifstream sersicdf("dfsersic.dat");
    
    for (int i = 0; i < n_psi; ++i)
    {
        denspsibulge >> Table_E[i] >> Dens_Psi_Bulge[i];
        sersicdf >> Table_E[i] >> DF_Sersic[i];
    }
    
    denspsibulge.close();
    sersicdf.close();
}
        
void ReadDensPsiHalo(void)
{
    ifstream denspsihalo("denspsihalo.dat");
    ifstream halodf("dfhalo.dat");
    
    for (int i = 0; i < n_psi; ++i)
    {
        denspsihalo >> Table_E[i] >> Dens_Psi_Halo[i];
        halodf >> Table_E[i] >> DF_NFW[i];
    }
    
    denspsihalo.close();
    halodf.close();
}

void WriteBDat(void)
{
    ofstream harmonicsfile("b.dat", ios::out);
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << setprecision(12) << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Bulge_Dens[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }  
    }   

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Bulge_Pot[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Bulge_FR[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }
    
    harmonicsfile.close();
}

void WriteHDat(void)
{
    ofstream harmonicsfile("h.dat", ios::out);
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << setprecision(12) << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Halo_Dens[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }  
    }   

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Halo_Pot[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << Halo_FR[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }
    
    harmonicsfile.close();
}
 
void WriteDBHDat(void)
{
    ofstream massfile("mass.dat", ios::out);
    ofstream harmonicsfile("dbh.dat", ios::out);
             
    massfile << total_disk_mass << " " << disk_edge << "\n" << bulge_mass << " " 
             << bulge_edge << "\n" << halo_mass << " " << halo_edge << endl;

//     harmonicsfile << G.c_halo << " " << G.v_halo << " " << G.a_halo << " " 
//                   << G.n_sersic << " " << G.v_bulge << " " << G.a_bulge << " " 
//                   << dr << " " << nr << " " << l_max << endl;
//     harmonicsfile << psi_0 << " " << halo_const << " " << bulge_const << endl;
//     harmonicsfile << G.M_Disk[0] << " " << G.R_Disk[0] << " " << G.Out_Disk[0] << " " 
//                   << G.Z_Disk[0] << " " << G.Dr_Trunc[0] << endl;
//     harmonicsfile << psi_crit << " " << psi_0 - psi_d << " " << G.bh_mass << endl;
//     harmonicsfile << disk_flag << " " << bulge_flag << " " << halo_flag << " " 
//                   << smbh_flag << endl;

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << setprecision(12) << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << A_Dens[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }  
    }   

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << A_Pot[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << F_R[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }
    }

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        harmonicsfile << rad << "   ";

        for (int l = 0; l < l_max+1; l+=2)
        {
            harmonicsfile << FR2[l/2][i] << "   ";
            if (l == l_max)
                harmonicsfile << endl;
        }     
    }

    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        double psi = Pot(rad, 0);
        harmonicsfile << rad << "   " << DiskDens(rad, 0, psi) << endl;
    }

    cout << "Final model written to file dbh.dat" << endl;
    
    massfile.close();
    harmonicsfile.close();
}

// Original code pretends it's all A_Pot so that the pot function can be used
// Here we use all the different vector<vector>s defined in galaxy.h

void ReadFreqs(void)
{
    ifstream dbhdatin("dbh.dat");
	float junk;
    
    // KK: READ HARMONICS IN FILE FILE INTO APOT ARRAY
    // THESE CAN THEN BE USED BY POT(R,Z) TO GIVE THE CORRESPONDING FUNCTION
    for (int i = 0; i < nr; ++i)
    {
        dbhdatin >> junk;
        for (int l = 0; l < l_max+1; l+=2)
        {
            dbhdatin >> A_Dens[l/2][i];
        }
    }
    
    // READ IN HARMONICS OF THE TOTAL POTENTIAL, TABULATE POT ON MAJ & MIN AXES
    for (int i = 0; i < nr; ++i)
    {
        dbhdatin >> junk;
        for (int l = 0; l < l_max+1; l+=2)
        {
            dbhdatin >> A_Pot[l/2][i];
        }
    }
    
    // KK: READ IN HARMONICS OF RADIAL FORCE, TABULATE ROTATION CURVE
    for (int i = 0; i < nr; ++i)
    {
        dbhdatin >> junk;
        for (int l = 0; l < l_max+1; l+=2)
        {
            dbhdatin >> F_R[l/2][i];
        }
    }
    
    // KK: READ IN HARMONICS OF 2ND DERIVATIVE OF PSI, TABULATE IN PSI2
    for (int i = 0; i < nr; ++i)
    {
        dbhdatin >> junk;
        for (int l = 0; l < l_max+1; l+=2)
        {
            dbhdatin >> FR2[l/2][i];
        }
    }
    
    for (int i = 0; i < nr; ++i)
    {
        double rad = Radius[i];
        double pot = Pot(rad, 0);
        
        Pot_Major_Tot[i] = pot;
        Pot_Minor[i] = Pot(0, rad);
        VC_Tot[i] = sqrt(max(0.0, rad*PotVector(rad, 0, F_R)));
        Psi2[i] = PotVector(rad, 0, FR2);
    }
    
    // KK: READ IN SURFACE DENSITY OF THE DISK, TABULATE
    for (int i = 0; i < nr; ++i)
    {
        dbhdatin >> junk >> Surf_Den[i];
    }
    
    dbhdatin.close();
    
    if (halo_flag)
    {
        ifstream hdatin("h.dat");

        // KK: READ HALO POTENTIAL, TABULATE IT ON MAJOR & INCLINED AXES
        for (int i = 0; i < nr; ++i)
        {
            hdatin >> junk;
            for (int l = 0; l < l_max+1; l+=2)
            {
                hdatin >> Halo_Dens[l/2][i];
            }
        }

        for (int i = 0; i < nr; ++i)
        {
            hdatin >> junk;
            for (int l = 0; l < l_max+1; l+=2)
            {
                hdatin >> Halo_Pot[l/2][i];
            }
        }

        // KK:  READ HALO RADIAL FORCE, TABULATE ITS ROTATION CURVE AND VERT FREQ
        for (int i = 0; i < nr; ++i)
        {
            hdatin >> junk;
            for (int l = 0; l < l_max+1; l+=2)
            {
                hdatin >> Halo_FR[l/2][i];
            }
        }

        for (int i = 0; i < nr; ++i)
        {
            double rad = Radius[i];
            double pot = PotVector(rad, 0, Halo_Pot);

            Pot_Major[i] = pot;
            Pot_Up[i] = PotVector(rad*cos(0.05), rad*sin(0.05), Halo_Pot);
            VC_Major[i] = sqrt(max(0.0, rad*PotVector(rad, 0, Halo_FR)));
            if (i != 0)
            {
                Vert_Freq[i] = sqrt(max(0.0, VC_Major[i]*VC_Major[i]+
                               800*(Pot_Up[i]-Pot_Major[i])))/rad;
            }
            else
            {
                Vert_Freq[i] = 0;//cout<<Vert_Freq[i]<<endl;
            }
        }

        hdatin.close();
    }
    
    if (bulge_flag)
    {
        ifstream bdatin("b.dat");

        // KK: READ BULGE POTENTIAL, TABULATE
        for (int i = 0; i < nr; ++i)
        {
            bdatin >> Radius[i];
            for (int l = 0; l < l_max+1; l+=2)
            {
                bdatin >> Bulge_Dens[l/2][i];
            }
        }

        for (int i = 0; i < nr; ++i)
        {
            bdatin >> junk;
            for (int l = 0; l < l_max+1; l+=2)
            {
                bdatin >> Bulge_Pot[l/2][i];
            }
        }

        // KK: READ BULGE RADIAL FORCE, MAKE ROT CURVE AND VERT FREQS. 
        for (int i = 0; i < nr; ++i)
        {
            bdatin >> junk;
            for (int l = 0; l < l_max+1; l+=2)
            {
                bdatin >> Bulge_FR[l/2][i];
            }
        }

        for (int i = 0; i < nr; ++i)
        {
            double rad = Radius[i];
            double pot = PotVector(rad, 0, Bulge_Pot);

            Pot_Major[i] = pot;
            Pot_Up[i] = PotVector(rad*cos(0.05), rad*sin(0.05), Bulge_Pot);
            VC_Bulge[i] = sqrt(max(0.0, rad*PotVector(rad, 0, Bulge_FR)));
            if (i != 0)
            {
                Vert_Freq_Bulge[i] = sqrt(max(0.0, VC_Bulge[i]*VC_Bulge[i]+
                                     800*(Pot_Up[i]-Pot_Major[i])))/rad;
            }
            else
            {
                Vert_Freq_Bulge[i] = 0;//cout<<Vert_Freq[i]<<endl;
            }
        }

        bdatin.close();
    }
    
// Output the data to freqdbh.dat. The columns are radius, Omega_H, Nu_H,
// Sigma_D, VC_tot, VC_b, Nu_B, PsiMaj_Tot, Psi
    
    ofstream freqout("freqdbh.dat");
    
    freqout << "0 " << VC_Major[0]/dr << " " << (4*Vert_Freq[1]-Vert_Freq[2])/3 << " "
            << Surf_Den[0] << " " << VC_Tot[0] << " " << VC_Bulge[0] << " "
            << (4*Vert_Freq_Bulge[1]-Vert_Freq_Bulge[2])/3 << " " << Pot_Major_Tot[0]
            << " " << Psi2[0] << endl;
    
    for (int i = 1; i < nr; ++i)
    {
        double rad = Radius[i];
        freqout << setw(8) << rad << " " << VC_Major[i]/rad << " " << Vert_Freq[i] << " "
                << Surf_Den[i] << " " << VC_Tot[i] << " " << VC_Bulge[i] << " "
                << Vert_Freq_Bulge[i] << " " << Pot_Major_Tot[i] << " " << Psi2[i] << endl;
    }
    
    freqout.close();
}
        
    
