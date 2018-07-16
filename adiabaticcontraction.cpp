//Contract the halo adiabatically. Note that the 'virial radius' of the halo
//is the halo radius from in.dbh plus 2*truncation_width. This should be
//treated as an outer radius beyond which baryon condensation does not occur.

#include "galaxy.h"

void ContractHalo(void)
{
    double max_halo_radius = r_max;//G.c_halo+2*G.drtrunc_halo;//r_max;
    double dr_ac = max_halo_radius/(nr_ac-1);
    double real_total_mass = 0;
    
    double Halo_Original_Dens[nr_ac];
    double Halo_Original_Mass[nr_ac];
    
    ofstream acfile, acfile2, fout;
    
    if (do_file_io)
    {
        acfile.open("adiabaticcontraction.dat", ios::out);
        acfile2.open("adiabaticcontraction.dat2", ios::out);
        fout.open("fout.dat", ios::out);
    }
    
    for(int i = 0; i < nr_ac; ++i)
    {
        Halo_Original_Radius[i] = i*dr_ac;//pow(10,rad);
    }
    
    for (int i = 1; i < nr_ac; ++i)
    {
        double radius2 = Halo_Original_Radius[i];
        double radius1 = Halo_Original_Radius[i-1];
        
        double rho = RawHaloProfile(radius2);//rho0 * exp(-2.0/alpha * (pow(radius/r_h, alpha)-1.0));
        double rho_p = RawHaloProfile(radius1);//rho0 * exp(-2.0/alpha * (pow((radius-frac*r_h)/r_h, alpha)-1.0));

        double volume = fourpi/3.0*(pow(radius2,3) - pow(radius1,3)); 

        //Use the midpoint density instead of the outer boundary density of each shell
        if (i > 1) 
            rho = (rho+rho_p)*0.5;
        
        Halo_Original_Dens[i] = rho;

        real_total_mass += rho*volume;

        Halo_Original_Mass[i] = real_total_mass;
        
//         double average_rho = total_mass_i / (fourpi/3.0*(pow(radius,3)));
//         if (average_rho < 200*rho_crit)
//         {
//             real_virial_radius = radius;
//             if (virial_radius > radius) virial_radius = radius;
//             break;
//         }
//         else if (radius > 998*r_h)
//         {
//             real_virial_radius = radius;
//             if (virial_radius > radius) virial_radius = radius;
//             break;
//         }  
    }
    
    if (disk_flag)
    {
        DiskMassWithRadius();
        DiskMassRPrime();
        //real_total_mass += Disk_Mass_Radius[nr_ac-1];
    }
    
    if (bulge_flag)
    {
        BulgeMassWithRadius();
        BulgeMassRPrime();
        cout << "rho0 " << G.rho_0 << " " << G.b_n << endl;
        //real_total_mass += Bulge_Mass_Radius[nr_ac-1];
    }
    
    baryon_frac = (Bulge_Mass_Radius[nr_ac-1]+Disk_Mass_Radius[nr_ac-1])/
                  real_total_mass;
    
    if (baryon_frac > 1)
    {
        cout << "Baryon mass is larger than halo mass! Lower M_disk, v_bulge, or increase v_halo\n"
             << Bulge_Mass_Radius[nr_ac-1] << " " << Disk_Mass_Radius[nr_ac-1] << " "
             << real_total_mass << endl;
        //cout << "Exiting..." << endl;
        baryon_frac = 1;
        exit(1);
    }
    else
    {
        cout << "Baryon fraction = " << baryon_frac << endl;
        cout << Bulge_Mass_Radius[nr_ac-1] << " " << Disk_Mass_Radius[nr_ac-1] 
             << " " << real_total_mass << " " << nr_ac << endl;
        cout << "Now getting contracted halo profile" << endl;
    }
    
    for (int i = 1; i < nr_ac; ++i)
    {
        acfile2 << Halo_Original_Radius[i]/r_max << " " << Halo_Original_Mass[i]*baryon_frac/real_total_mass
                << " " << Halo_Original_Mass[i]*(1-baryon_frac)/real_total_mass << " "
                << Halo_Original_Radius[i]/r_max << " " << (Bulge_Mass_Radius[i] + Disk_Mass_Radius[i])/real_total_mass
                << " 0.0" << endl;
//         acfile2 << Halo_Original_Radius[i]/r_max << " " << Halo_Original_Mass[i]*baryon_frac
//                 << " " << Halo_Original_Mass[i]*(1-baryon_frac) << " "
//                 << Halo_Original_Radius[i]/r_max << " " << Bulge_Mass_Radius[i] + Disk_Mass_Radius[i]
//                 << " 0.0" << endl;
    }
                   
    //double rhalo_mass = real_total_mass*(1-baryon_frac);
    double disk_mass, disk_mass_prime, bulge_mass, bulge_mass_prime;
    double mass_i = 0, mass_i_org = 0;
    double r_1, r_1_old;
    
    //ofstream fout("fout.dat", ios::out);
    if (contraction_prescription == 1)
    {
        for (int i = 1; i < nr_ac; ++i)
        {
            Halo_Original_Dens[i] *= (1-baryon_frac);
            
            double radius2 = Halo_Original_Radius[i], radius1;
            
	        if (i==0) 
            {
                radius1 = 0;
            }
            else
            {
                radius1 = Halo_Original_Radius[i-1];
            }
	        //if (i==0) radius1 = 0;

            double rho = RawHaloProfile(radius2);
            double rho_p = RawHaloProfile(radius1);
	        double volume = fourpi/3.0 * (pow(radius2,3) - pow(radius1,3)); 
            
            if (i > 1)
                rho = (rho+rho_p)*0.5;
            
	        mass_i += rho*volume;
    
	        double invariant = radius2*mass_i;
	        
	        //get m_f for the final halo configuration
	        //we don't know the radius yet - plug this into the adiabatic 
	        //contraction equation to get the radius
	        //Halo_Mass_Radius[(int)(radius/frac/r_h+0.01)] = mass_i * (1.0-baryon_frac);
	        Halo_Mass_Radius[i] = mass_i*(1.0-baryon_frac);
	        
	        if (i == 1)
	        {
	            r_1_old = 0;
	            Halo_Mass_Radius[0] = 0;
	        }
	        else
	        {
	            r_1_old = r_1;
	        }
            
	        //Newton-Raphson method to get the radius
	        //the function is the adiabatic contraction expression set to zero
	        double r_0 = 0.5 * radius2; //guess the first value
	        r_1 = 0.5 * radius2;
	        
	        //Get sigma (which is really M_disk_tot) by evaluating the mass inside
	        //the virial radius and multiplying by F
	        //double sigma = G.M_Disk[j];//total_mass_i*(baryon_frac);
	        double f_of_r, fprime_of_r;
	        int nrcount = 0;
	        
	        do
	        {
	            r_0 = r_1;//double sigma = G.M_Disk[0], rscale = G.R_Disk[0], radius = radius2;
                bulge_mass = GetBulgeMass(r_0);
                bulge_mass_prime = GetBulgeMassPrime(r_0);
                disk_mass = GetDiskMassR(r_0);
                disk_mass_prime = GetDiskMassPrime(r_0);
                f_of_r = r_0*(1-baryon_frac)*mass_i + r_0*disk_mass +
                         r_0*bulge_mass - radius2*mass_i;
                fprime_of_r = (1-baryon_frac)*mass_i + disk_mass + r_0*disk_mass_prime +
                              bulge_mass + r_0*bulge_mass_prime;
//                 f_of_r = r_0*sigma + r_0*(1-baryon_frac)*mass_i - r_0*sigma*exp(-r_0/rscale) - 
//                          r_0*r_0/rscale*sigma*exp(-r_0/rscale) - radius*mass_i;// + r_0*bulge_mass;
//                 fprime_of_r = sigma + (1-baryon_frac)*mass_i - sigma*exp(-r_0/rscale) -
//                               r_0/rscale*sigma*exp(-r_0/rscale) +
//                               r_0*r_0/rscale/rscale*sigma*exp(-r_0/rscale);// +
//                               //bulge_mass + r_0*bulge_mass_prime;
	            r_1 = r_0 - f_of_r / fprime_of_r;
	            ++nrcount;
	            //cout << r_1<<" " << r_0<<" "<<nrcount<<endl;
	        }
	        while(fabs(r_1 - r_0) > 0.000000001 && nrcount < 100);
	         
	        double rho_f = (Halo_Mass_Radius[i] - Halo_Mass_Radius[i-1]) / 
                           (fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3)));
            
            double density_conversion = pow(r_max,3)/real_total_mass;
            
            if (do_file_io)
            {
//                 acfile << radius2 << " " << r_1 << "    " << Halo_Mass_Radius[i]
//                        << "    " << rho << " " << rho_f << "     " << nrcount 
//                        << "    " << log10(radius2) << " " << log10(r_1) << " "
//                        << log10(rho) << " " << log10(rho_f) << "        " 
//                        << Disk_Mass_Radius[i] << " " << invariant << " " 
//                        << r_1*Halo_Mass_Radius[i] + r_1*GetDiskMassR(r_1) + r_1*GetBulgeMass(r_1) 
//                        << " " << rho_f*(fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3))) << endl;

                acfile << radius2/r_max << " " << r_1/r_max << "    " << Halo_Mass_Radius[i]
                       << "    " << rho*density_conversion << " " << rho_f*density_conversion << "     " << nrcount 
                       << "    " << log10(radius2) << " " << log10(r_1) << " "
                       << log10(rho) << " " << log10(rho_f) << "        " 
                       << Disk_Mass_Radius[i] << " " << invariant << " " 
                       << r_1*Halo_Mass_Radius[i] + r_1*GetDiskMassR(r_1) + r_1*GetBulgeMass(r_1) 
                       << " " << radius2/r_1 << " " << rho_f/rho << endl;
            }
            
            Halo_AC_Radius[i] = r_1;
            Halo_AC_Dens[i] = rho_f;
            
            //cout << "Done step " << i << endl;
        }
    }
	else if (contraction_prescription == 2)
    {
	    double A = 0.85, w = 0.8;
	    A *= pow(max_halo_radius, 1-w);    
	        
	    for (int i = 1; i < nr_ac; ++i)
	    {
            double radius2 = Halo_Original_Radius[i];
            double radius1 = Halo_Original_Radius[i-1];
	        if (i==0) radius1 = 0;
	        
	        double radius_oa2 = A*pow(radius2, w);
	        double radius_oa1 = A*pow(radius1, w);
	       
	        double rho = RawHaloProfile(radius_oa2);
	        double rho_p = RawHaloProfile(radius_oa1);
	       
	        double volume = fourpi/3.0*(pow(radius_oa2,3) - pow(radius_oa1,3)); 
            
            if (i > 1)
                rho = (rho+rho_p)*0.5;
            
	        mass_i += rho*volume;
	        //invariant = radius*mass_i;
	         
	        double rho_org = RawHaloProfile(radius2);
	        double rho_org_p = RawHaloProfile(radius1);
	       
	        double volume_org = fourpi/3.0*(pow(radius2,3) - pow(radius1,3)); 
            
            if (i > 1)
                rho_org = (rho_org+rho_org_p)*0.5;
            
	        mass_i_org += rho_org*volume_org;
	        //invariant = radius*mass_i;
	       
	        //get m_f for the final halo configuration
	        //we don't know the radius yet - plug this into the adiabatic 
	        //contraction equation to get the radius
	        Halo_Mass_Radius[i] = mass_i_org * (1.0-baryon_frac);
	        double mass_ii = mass_i * (1.0-baryon_frac);
	       
	        if (i == 1)
	        {
	            r_1_old = 0;
	            Halo_Mass_Radius[0] = 0;
	        }
	        else
	        {
	            r_1_old = r_1;
	        }
	        
	        //Newton-Raphson method to get the radius
	        //the function is the adiabatic contraction expression set to zero
	        double r_0 = 0.5*radius2; //guess the first value
	        r_1 = 0.5*radius2;
	        
	        //Get sigma (which is really M_disk_tot) by evaluating the mass inside
	        //the virial radius and multiplying by F
	        //double sigma = G.M_Disk[j];//total_mass_i*(baryon_frac);
	        double f_of_r, fprime_of_r;
	        int nrcount = 0;
	        
	        do
	        {
	            r_0 = r_1;
                double r_0_oa = A*pow(r_0, w);//if (r_0<0) r_0_oa = 0;
                double drdr = w*r_0_oa/r_0;//cout << "dr " << drdr << endl;
                bulge_mass = GetBulgeMass(r_0_oa);//cout << "b " << bulge_mass << endl;
                bulge_mass_prime = GetBulgeMassPrime(r_0_oa);//cout << "bb " << bulge_mass_prime << endl;
                disk_mass = GetDiskMassR(r_0_oa);//cout << "d " << disk_mass << endl;
                disk_mass_prime = GetDiskMassPrime(r_0_oa);//cout << "dd " << disk_mass_prime << endl;
                f_of_r = r_0*(1-baryon_frac)*mass_i + r_0*disk_mass + 
                         r_0*bulge_mass - radius2*mass_i;
                fprime_of_r = (1-baryon_frac)*mass_i + disk_mass + r_0*disk_mass_prime*drdr +
                              bulge_mass + r_0*bulge_mass_prime*drdr;
// 	            f_of_r = r_0*sigma + r_0*(1-baryon_frac)*mass_i - r_0*sigma*exp(-A*pow(r_0, w)/rscale) - 
// 	                     A*pow(r_0, w+1)/rscale*sigma*exp(-A*pow(r_0, w)/rscale) - radius*mass_i;
// 	            fprime_of_r = sigma + (1-baryon_frac)*mass_i - sigma*exp(-A*pow(r_0, w)/rscale) -
// 	                          A*pow(r_0, w)/rscale*sigma*exp(-A*pow(r_0, w)/rscale) +
// 	                          A*A*w*pow(r_0, 2*w)/rscale/rscale*sigma*exp(-A*pow(r_0, w)/rscale);
	            r_1 = r_0 - f_of_r / fprime_of_r;
	            ++nrcount;
	            //cout << r_1 << " " << r_0 << " " << r_0_oa << " " << nrcount<<endl;
	        }
	        while(fabs(r_1 - r_0) > 0.0000000001 && nrcount < 100);
	            
	        double rho_f = (Halo_Mass_Radius[i] - Halo_Mass_Radius[i-1]) / 
	                       (fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3)));
	                
	        //radius_oa = A*pow(r_1, w);
            
            double density_conversion = 1;//pow(r_max,3)/rhalo_mass;
            
            if (do_file_io)
            {
//                 acfile << radius2 << " " << r_1 << "    " << Halo_Mass_Radius[i]
//                        << "    " << rho << " " << rho_f << "     " << nrcount 
//                        << "    " << log10(radius2) << " " << log10(r_1) << " "
//                        << log10(rho) << " " << log10(rho_f) << "        " 
//                        << Disk_Mass_Radius[i] << " " << invariant << " " 
//                        << r_1*Halo_Mass_Radius[i] + r_1*GetDiskMassR(r_1) + r_1*GetBulgeMass(r_1) 
//                        << " " << rho_f*(fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3))) << endl;

                acfile << radius2/r_max << " " << r_1/r_max << "    " << Halo_Mass_Radius[i]
                       << "    " << rho*density_conversion << " " << rho_f*density_conversion << "     " << nrcount 
                       << "    " << log10(radius2) << " " << log10(r_1) << " "
                       << log10(rho) << " " << log10(rho_f) << "        " 
                       << Disk_Mass_Radius[i] << " " << radius2/r_1 << " " << rho_f/rho << endl;
            }
            
            Halo_AC_Radius[i] = r_1;
            Halo_AC_Dens[i] = rho_f;
	    }
    }
	else if (contraction_prescription == 3)
	{
	    for (int i = 1; i < nr_ac; i++)
	    {
            double radius2 = Halo_Original_Radius[i];
            double radius1 = Halo_Original_Radius[i-1];
	        if (i==0) radius1 = 0;
	        
	        double rho =  RawHaloProfile(radius2);
	        
	        double volume = fourpi/3.0*(pow(radius2,3) - pow(radius1,3)); 
	        mass_i += rho*volume;
	        //invariant = radius*mass_i;
	        
	        //get m_f for the final halo configuration
	        //we don't know the radius yet - plug this into the adiabatic 
	        //contraction equation to get the radius
	        Halo_Mass_Radius[i] = mass_i*(1.0-baryon_frac);
	        
	        if (i == 1)
	        {
	            r_1_old = 0;
	            Halo_Mass_Radius[0] = 0;
	        }
	        else
	        {
	            r_1_old = r_1;
	        }
	        
	        //Newton-Raphson method to get the radius
	        //the function is the adiabatic contraction expression set to zero
	        double r_0 = 0.5*radius2; //guess the first value
	        r_1 = 0.5*radius2;
	        
	        //Get sigma (which is really M_disk_tot) by evaluating the mass inside
	        //the virial radius and multiplying by F
	        //double sigma = total_mass_i*(baryon_frac);
	        double f_of_r, fprime_of_r;
	        int nrcount = 0;

	        do
	        {
	            r_0 = r_1;
                bulge_mass = GetBulgeMass(r_0);
                bulge_mass_prime = GetBulgeMassPrime(r_0);
                disk_mass = GetDiskMassR(r_0);
                disk_mass_prime = GetDiskMassPrime(r_0);
                f_of_r = r_0*(1-baryon_frac)*mass_i + r_0*disk_mass + 
                         r_0*bulge_mass - radius2*mass_i;
                fprime_of_r = (1-baryon_frac)*mass_i + disk_mass + r_0*disk_mass_prime +
                              bulge_mass + r_0*bulge_mass_prime;
// 	            f_of_r = r_0*sigma + r_0*(1-baryon_frac)*mass_i - r_0*sigma*exp(-r_0/rscale) - 
// 	                     r_0*r_0/rscale*sigma*exp(-r_0/rscale) - radius*mass_i;
// 	            fprime_of_r = sigma + (1-baryon_frac)*mass_i - sigma*exp(-r_0/rscale) -
// 	                          r_0/rscale*sigma*exp(-r_0/rscale) +
// 	                          r_0*r_0/rscale/rscale*sigma*exp(-r_0/rscale);
	            r_1 = r_0 - f_of_r / fprime_of_r;
	            ++nrcount;
	            //cout << r_1<<" " << r_0<<" "<<nrcount<<endl;
	        }
	        while(fabs(r_1 - r_0) > 0.000001);
	            
	        double r_2, r_2_old;
	        
	        if (i == 0)
	        {
	            r_2_old = 0;
	        }
	        else
	        {
	            r_2_old = r_2;
	        }
	        
	        double a = 0.3, n = 2;
	        
	        double massratio = mass_i/(Halo_Mass_Radius[i] + GetDiskMassR(r_1) +
                               GetBulgeMass(r_1));
	        
	        double radratio = 1 + a*(pow(massratio, n) - 1);
	        r_2 = radratio*radius2;
	
	        double rho_f = (Halo_Mass_Radius[i] - Halo_Mass_Radius[i-1]) / 
	                       (fourpi/3.0*(pow(r_2,3) - pow(r_2_old,3)));
	        
            double density_conversion = 1;//pow(r_max,3)/rhalo_mass;
            
            if (do_file_io)
            {
//                 acfile << radius2 << " " << r_1 << "    " << Halo_Mass_Radius[i]
//                        << "    " << rho << " " << rho_f << "     " << nrcount 
//                        << "    " << log10(radius2) << " " << log10(r_1) << " "
//                        << log10(rho) << " " << log10(rho_f) << "        " 
//                        << Disk_Mass_Radius[i] << " " << invariant << " " 
//                        << r_1*Halo_Mass_Radius[i] + r_1*GetDiskMassR(r_1) + r_1*GetBulgeMass(r_1) 
//                        << " " << rho_f*(fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3))) << endl;

                acfile << radius2/r_max << " " << r_1/r_max << "    " << Halo_Mass_Radius[i]
                       << "    " << rho*density_conversion << " " << rho_f*density_conversion << "     " << nrcount 
                       << "    " << log10(radius2) << " " << log10(r_1) << " "
                       << log10(rho) << " " << log10(rho_f) << "        " 
                       << Disk_Mass_Radius[i] << endl;
            }
            
            Halo_AC_Radius[i] = r_1;
            Halo_AC_Dens[i] = rho_f;
	    }
    }
    if (contraction_prescription == 0)
    {
        for (int i = 1; i < nr_ac; ++i)
        {
            Halo_Original_Dens[i] *= (1-baryon_frac);
            
            double radius2 = Halo_Original_Radius[i], radius1;
	        if (i==0) 
            {
                radius1 = 0;
            }
            else
            {
                radius1 = Halo_Original_Radius[i-1];
            }
	        //if (i==0) radius1 = 0;

            double rho = RawHaloProfile(radius2);
            double rho_p = RawHaloProfile(radius1);
	        double volume = fourpi/3.0 * (pow(radius2,3) - pow(radius1,3)); 
            
            if (i > 1)
                rho = (rho+rho_p)*0.5;
            
	        mass_i += rho*volume;
            
	        Halo_Mass_Radius[i] = mass_i*(1.0-baryon_frac);

	        if (i == 1)
	        {
	            r_1_old = 0;
	            Halo_Mass_Radius[0] = 0;
	        }
	        else
	        {
	            r_1_old = r_1;
	        }
            
            r_1 = radius2;
            
            int nrcount = 1;
	        
	        double rho_f = (Halo_Mass_Radius[i] - Halo_Mass_Radius[i-1]) / 
                           (fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3)));
            
            double density_conversion = 1;//pow(r_max,3)/rhalo_mass;
            
            if (do_file_io)
            {
//                 acfile << radius2 << " " << r_1 << "    " << Halo_Mass_Radius[i]
//                        << "    " << rho << " " << rho_f << "     " << nrcount 
//                        << "    " << log10(radius2) << " " << log10(r_1) << " "
//                        << log10(rho) << " " << log10(rho_f) << "        " 
//                        << Disk_Mass_Radius[i] << " " << invariant << " " 
//                        << r_1*Halo_Mass_Radius[i] + r_1*GetDiskMassR(r_1) + r_1*GetBulgeMass(r_1) 
//                        << " " << rho_f*(fourpi/3.0*(pow(r_1,3) - pow(r_1_old,3))) << endl;

                acfile << radius2/r_max << " " << r_1/r_max << "    " << Halo_Mass_Radius[i]
                       << " " << Halo_Original_Dens[i]*density_conversion 
                       << "    " << rho*density_conversion << " " << rho_f*density_conversion << "     " << nrcount 
                       << "    " << log10(radius2) << " " << log10(r_1) << " "
                       << log10(rho) << " " << log10(rho_f) << "        " 
                       << Disk_Mass_Radius[i] << endl;
            }
            
            Halo_AC_Radius[i] = r_1;
            Halo_AC_Dens[i] = rho_f;
        }
    }
    
    //Now interpolate to get the density at each grid point 
    //(rather than the midpoint of each bin).
    for (int i = 1; i < nr_ac-1; ++i)
    {
        double dens = Halo_AC_Dens[i]+Halo_AC_Dens[i+1];
        Halo_AC_Dens[i] = dens*0.5;
    }
    
    cout << "Got contraction" << endl;
    
    Halo_AC_Dens[0] = 2*Halo_AC_Dens[1] - Halo_AC_Dens[2];
    Halo_AC_Dens_D[nr_ac-1] = 2*Halo_AC_Dens_D[nr_ac-2] - Halo_AC_Dens_D[nr_ac-3];
    
    //Get derivatives of the new halo density
    for (int i = 1; i < nr_ac-1; ++i)
    {
        double drho = Halo_AC_Dens[i+1]-Halo_AC_Dens[i-1];
        double dr1 = Halo_AC_Radius[i+1]-Halo_AC_Radius[i-1];
        
        Halo_AC_Dens_D[i] = drho/dr1;
    }
    
    Halo_AC_Dens_D[0] = 2*Halo_AC_Dens_D[1] - Halo_AC_Dens_D[2];
    Halo_AC_Dens_D[nr_ac-1] = 2*Halo_AC_Dens_D[nr_ac-2] - Halo_AC_Dens_D[nr_ac-3];
    
    for (int i = 1; i < nr_ac-1; ++i)
    {
        double drho = Halo_AC_Dens_D[i+1]-Halo_AC_Dens_D[i-1];
        double dr1 = Halo_AC_Radius[i+1]-Halo_AC_Radius[i-1];
        
        Halo_AC_Dens_DD[i] = drho/dr1;
    }
    
    Halo_AC_Dens_DD[0] = 2*Halo_AC_Dens_DD[1] - Halo_AC_Dens_DD[2];
    Halo_AC_Dens_DD[nr_ac-1] = 2*Halo_AC_Dens_DD[nr_ac-2] - Halo_AC_Dens_DD[nr_ac-3];
    
    if (do_file_io)
    {
//        for (int i = 0; i < nr; ++i)
//         {
//             //acfile2 << i << " " << Radius[i] << " " << log10(Radius[i]) << " "
//             //        << log10(HaloProfileDens(Radius[i])) << endl;
//             acfile2 << Radius[i] << " " << setprecision(12) << HaloProfileDens(Radius[i]) << " "
//                     << HaloProfileDensPrime(Radius[i]) << " " 
//                     << HaloProfileDens2Prime(Radius[i]) << " " << endl;
//         }
//         
            
        double density_conversion = pow(r_max,3)/real_total_mass;
        
        for (int i = 1; i < nr_ac; ++i)
        {
            double rad = Halo_AC_Radius[i], dens = Halo_AC_Dens[i];
            double dens2 = HaloProfileDens(rad);
            
//             fout << Halo_Original_Radius[i]/r_max << " " << rad/r_max << " " 
//                  << dens << " " << dens2 << "      " 
//                  << log10(rad/r_max) << " " << log10(dens*rad*rad/r_max/r_max) << " "
//                  << Halo_AC_Dens_D[i] << " " << HaloProfileDensPrime(rad) << "    "
//                  << Halo_AC_Dens_DD[i] << " " << HaloProfileDens2Prime(rad) << endl;
            fout << Halo_Original_Radius[i]/r_max << " " << rad/r_max << " " 
                 << dens*density_conversion << " " << dens2*density_conversion << "      " 
                 << log10(rad/r_max) << " " << log10(dens*rad*rad*density_conversion/r_max/r_max) << " "
                 << Halo_AC_Dens_D[i] << " " << HaloProfileDensPrime(rad) << "    "
                 << Halo_AC_Dens_DD[i] << " " << HaloProfileDens2Prime(rad) << "    "
                 << endl;
        }
            
            //double rad = Radius[i];
            //fout << rad << " " << rad*(1-baryon_frac)*Halo_Mass_Radius[i]+rad*GetDiskMassR(rad)+
            //        rad*GetBulgeMass(rad) - Halo_Original_Radius[0]*Halo_Mass_Radius[0]/(1.0-baryon_frac) 
            //     << endl;
    }
    
}

//Return the interpolated derivative of the bulge mass function at an 
//arbitrary radius

double GetBulgeMassPrime(double r)
{
    if (!bulge_flag) 
    {
        return 0;
    }
    
    if (r < 0) 
    {
        return -gsl_spline_eval(bulge_mass_prime_spline, -r, bulge_mass_prime_acc);
    }
    else if (r > Halo_Original_Radius[nr_ac-1])
    {
        return gsl_spline_eval(bulge_mass_prime_spline, Halo_Original_Radius[nr_ac-1], bulge_mass_prime_acc);
    }
    else
    {
        return gsl_spline_eval(bulge_mass_prime_spline, r, bulge_mass_prime_acc);
    }
    
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Getbulgerprime finds out of range indices. Exiting..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        exit(1);
    }
    else if (ihi > nr_ac-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_ac-1;
    }
    
    double r1 = Halo_Original_Radius[ihi-1];
    double r2 = Halo_Original_Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*Bulge_Mass_Prime[ihi] + tm1*Bulge_Mass_Prime[ihi-1];
}    
    
//Return the interpolated bulge mass at an arbitrary radius

double GetBulgeMass(double r)
{
    if (!bulge_flag) 
    {
        return 0;
    }
    
    if (r < 0) 
    {
        return -gsl_spline_eval(bulge_mass_spline, -r, bulge_mass_acc);
    }
    else if (r > Halo_Original_Radius[nr_ac-1])
    {
        return gsl_spline_eval(bulge_mass_spline, Halo_Original_Radius[nr_ac-1], bulge_mass_acc);
    }
    else
    {
        return gsl_spline_eval(bulge_mass_spline, r, bulge_mass_acc);
    }
    
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Getbulgerprime finds out of range indices. Exiting..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        exit(1);
    }
    else if (ihi > nr_ac-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_ac-1;
    }
    
    double r1 = Halo_Original_Radius[ihi-1];
    double r2 = Halo_Original_Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*Bulge_Mass_Radius[ihi] + tm1*Bulge_Mass_Radius[ihi-1];
}
    
//Get the bulge mass derivative as a function of radius, then initialise a 
//spline to interpolate

void BulgeMassRPrime(void)
{
    for (int i = 1; i < nr_ac-1; ++i)
    {
        double dmass = Bulge_Mass_Radius[i+1]-Bulge_Mass_Radius[i-1];
        double dr1 = Halo_Original_Radius[i+1]-Halo_Original_Radius[i-1];
        
        Bulge_Mass_Prime[i] = dmass/dr1;
        //cout << i << " " << Bulge_Mass_Prime[i] << endl;
    }
    
    Bulge_Mass_Prime[0] = max(0.0, 2*Bulge_Mass_Prime[1] - Bulge_Mass_Prime[2]);
                          //(Bulge_Mass_Radius[1]-Bulge_Mass_Radius[0])/
                          //(Halo_Original_Radius[1]-Halo_Original_Radius[0]);
    Bulge_Mass_Prime[nr_ac-1] = 2*Bulge_Mass_Prime[nr_ac-2] - Bulge_Mass_Prime[nr_ac-3];
    
    bulge_mass_prime_acc = gsl_interp_accel_alloc();
    bulge_mass_prime_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(bulge_mass_prime_spline, &Halo_Original_Radius[0], &Bulge_Mass_Prime[0], nr_ac);
}

//Get the bulge mass as a function of radius, then initialise a spline
//to interpolate

void BulgeMassWithRadius(void)
{
    double max_bulge_radius = r_max;//50*G.a_bulge;
    double bulge_mass = 0;
    
    Bulge_Mass_Radius.at(0) = 0;
    
    for (int i = 1; i < nr_ac; ++i)
    {
        double radius2 = Halo_Original_Radius[i];
        double radius1 = Halo_Original_Radius[i-1];
        
        double rho = SersicDens(radius2);
        double rho_p = SersicDens(radius1);

        double volume = fourpi/3.0 * (pow(radius2,3) - pow(radius1,3)); 

        //Use the midpoint density instead of the outer boundary density of each shell
        if (i > 1) 
            rho = (rho+rho_p)*0.5;

        bulge_mass += rho*volume;
        
        Bulge_Mass_Radius[i] = bulge_mass;
        //cout << i << " " << bulge_mass << endl;
    }
    
    bulge_mass_acc = gsl_interp_accel_alloc();
    bulge_mass_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(bulge_mass_spline, &Halo_Original_Radius[0], &Bulge_Mass_Radius[0], nr_ac);
}

//Return the interpolated derivative of the disk mass function at an 
//arbitrary radius

double GetDiskMassPrime(double r)
{
    if (!disk_flag) 
    {
        return 0;
    }
    
    if (r < 0) 
    {
        return -gsl_spline_eval(disk_mass_prime_spline, -r, disk_mass_prime_acc);
    }
    else if (r > Halo_Original_Radius[nr_ac-1])
    {
        return gsl_spline_eval(disk_mass_prime_spline, Halo_Original_Radius[nr_ac-1], disk_mass_prime_acc);
    }
    else
    {
        return gsl_spline_eval(disk_mass_prime_spline, r, disk_mass_prime_acc);
    }
    
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Getdiskrprime finds out of range indices. Exiting..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        exit(1);
    }
    else if (ihi > nr_ac-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_ac-1;
    }
    
    double r1 = Halo_Original_Radius[ihi-1];
    double r2 = Halo_Original_Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*Disk_Mass_Prime[ihi] + tm1*Disk_Mass_Prime[ihi-1];
}    
    
//Return the interpolated disk mass at an arbitrary radius

double GetDiskMassR(double r)
{
    if (!disk_flag) 
    {
        return 0;
    }
    
    if (r < 0) 
    {
        return -gsl_spline_eval(disk_mass_spline, -r, disk_mass_acc);
    }
    else if (r > Halo_Original_Radius[nr_ac-1])
    {
        return gsl_spline_eval(disk_mass_spline, Halo_Original_Radius[nr_ac-1], disk_mass_acc);
    }
    else
    {
        return gsl_spline_eval(disk_mass_spline, r, disk_mass_acc);
    }
        
    int ihi = ceil(r/dr);
    
    if(r < dr)
    {
        ihi = 1;
    }
    else if (ihi < 1)
    {
        cout << "Getdiskrprime finds out of range indices. Exiting..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        exit(1);
    }
    else if (ihi > nr_ac-1)
    {
        //cout << "GettotalPsi finds out of range indices. Continuing..." << endl;
        //cout <<"ihi "<<ihi<<" "<<r<<" "<<log_r<<" "<<log_dr<<" "<<log_rmax << endl;
        ihi = nr_ac-1;
    }
    
    double r1 = Halo_Original_Radius[ihi-1];
    double r2 = Halo_Original_Radius[ihi];
    double t = (r-r1)/(r2-r1);
    double tm1 = 1-t;
    
    return t*Disk_Mass_Radius[ihi] + tm1*Disk_Mass_Radius[ihi-1];
}
    
void DiskMassRPrime(void)
{
    for (int i = 1; i < nr_ac-1; ++i)
    {
        double dmass = Disk_Mass_Radius[i+1]-Disk_Mass_Radius[i-1];
        double dr1 = Halo_Original_Radius[i+1]-Halo_Original_Radius[i-1];
        
        Disk_Mass_Prime[i] = dmass/dr1;
    }
    
    Disk_Mass_Prime[0] = 2*Disk_Mass_Prime[1] - Disk_Mass_Prime[2];
    Disk_Mass_Prime[nr_ac-1] = 2*Disk_Mass_Prime[nr_ac-2] - Disk_Mass_Prime[nr_ac-3];
    
    disk_mass_prime_acc = gsl_interp_accel_alloc();
    disk_mass_prime_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(disk_mass_prime_spline, &Halo_Original_Radius[0], &Disk_Mass_Prime[0], nr_ac);
}

//Get the disk mass as a function of radius, then initialise a spline
//to interpolate

void DiskMassWithRadius(void)
{
    double max_disk_radius = G.Out_Disk[0] + 2*G.Dr_Trunc[0];
    
    Disk_Mass_Radius.at(0) = 0;
    
    for (int i = 1; i < nr_ac; ++i)
    {
        double radius2 = Halo_Original_Radius[i];
        double radius1 = Halo_Original_Radius[i-1];
        double sigma = 0, sigma_p = 0;
        
        for (int j = 0; j < disk; ++j)
        {
            double trunc_fac1 = GetTrunc(radius1, j);
            double trunc_fac2 = GetTrunc(radius2, j);
		    double Sigma_Profile[2], Rho_Profile[2];
            double z = 0;

            DiskProfile(radius2, z, j, Sigma_Profile, Rho_Profile);
            
            sigma += Disk_Const[j]*Sigma_Profile[0]*trunc_fac2;

            DiskProfile(radius1, z, j, Sigma_Profile, Rho_Profile);
            
            sigma_p += Disk_Const[j]*Sigma_Profile[0]*trunc_fac1;
        }

        double area = PI*(pow(radius2,2) - pow(radius1,2)); 

        //Use the midpoint density instead of the outer boundary density of each shell
        if (i > 1) 
            sigma = (sigma+sigma_p)*0.5;

        disk_mass += sigma*area;
        
        Disk_Mass_Radius[i] = disk_mass;
    }
    
    disk_mass_acc = gsl_interp_accel_alloc();
    disk_mass_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(disk_mass_spline, &Halo_Original_Radius[0], &Disk_Mass_Radius[0], nr_ac);
}

            
void HaloMassDwithRadius(void)
{
    halo_dens_acc = gsl_interp_accel_alloc();
    halo_dens_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(halo_dens_spline, &Halo_AC_Radius[0], &Halo_AC_Dens[0], nr_ac);
    
    halo_dens_prime_acc = gsl_interp_accel_alloc();
    halo_dens_prime_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(halo_dens_prime_spline, &Halo_AC_Radius[0], &Halo_AC_Dens_D[0], nr_ac);
    
    halo_dens_2prime_acc = gsl_interp_accel_alloc();
    halo_dens_2prime_spline = gsl_spline_alloc(gsl_interp_cspline, nr_ac);
    
    gsl_spline_init(halo_dens_2prime_spline, &Halo_AC_Radius[0], &Halo_AC_Dens_DD[0], nr_ac);
}

            
