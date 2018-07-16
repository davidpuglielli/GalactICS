// Calculate the LOS quantities of a galaxy using the DF and compare to
// observations. Here we get circular velocity, LOS stellar velocity and
// dispersion, and surface brightness profile

// We need to read in the observational data earlier so as not to have to 
// do it at each step of a chain, as is done in the original code

#include "galaxy.h"
#include "chisq.h"

double GetChiSquare(vector<double> &Chi_Square)
{
    //Read in the astronomical parameters
    ifstream astro_in("in.astronomical");
    ifstream error_in("in.errorfactors");
    
    if (!(astro_in.is_open()&&error_in.is_open()))
    {
        cerr << "in.astronomical or in.errorfactors is missing. Exiting..." << endl;
    }
    
    Astro.ML_Disk.resize(disk);
    
    for (int i = 0; i < disk; ++i)
    {
        astro_in >> Astro.ML_Disk[i];
    }
    
    astro_in >> Astro.ml_bulge >> Astro.inclination 
             >> Astro.distance >> Astro.v_sys >> Astro.r_sys;
    
    Astro.inclination *= PI/180;
    
    cos_inclination = cos(Astro.inclination);
    sin_inclination = sin(Astro.inclination);
    
    for (int i = 0; i < Error_Factors.size(); ++i)
    {
        error_in >> Error_Factors[i];
    }
    
    for (int i = 0; i < DOF.size(); ++i)
    {
        DOF[i] = 0;
    }
    
    T = gsl_rng_mt19937;
    rand_gen = gsl_rng_alloc(T);
    gsl_rng_set(rand_gen, 0);
    
    //for (int i = 0; i < VC_Radii_Orig.size(); ++i)
    //{
    //    cout << VC_Radii_Orig[i] << endl;
    //}
    
    //Convert the observations to model units using the astronomical parameters
    ConvertRadius(VC_Radii_Orig, VC_Radii);
    ConvertRadius(SVel_Radii_Orig, SVel_Radii);
    ConvertRadius(SB_Radii_Orig, SB_Radii);
    ConvertRadius(SDisp_Radii_Orig, SDisp_Radii);
    //ConvertRadius(HIRot_Radii_Orig, HIRot_Radii);
    //ConvertRadius(HIDisp_Radii_Orig, HIDisp_Radii);
    
    ConvertVelocities(VC_Data_Orig, VC_Error_Orig, VC_Data, VC_Error);
    ConvertVelocities(SVel_Data_Orig, SVel_Error_Orig, SVel_Data, SVel_Error);
    ConvertVelocities(SDisp_Data_Orig, SDisp_Error_Orig, SDisp_Data, SDisp_Error);
    //ConvertVelocities(HIRot_Data_Orig, HIRot_Error_Orig, HIRot_Data, HIRot_Error);
    //ConvertVelocities(HIDisp_Data_Orig, HIDisp_Error_Orig, HIDisp_Data, HIDisp_Error);
    
    //Evaluate the chi squares
    Chi_Square[0] = GetSurfaceBrightness(SB_Radii, SB_Data, SB_Error);
    Chi_Square[1] = GetCircVelocity(VC_Radii, VC_Data, VC_Error);
    Chi_Square[2] = GetStarVelocity(SVel_Radii, SVel_Data, SVel_Error);
    Chi_Square[3] = GetStarDispersion(SDisp_Radii, SDisp_Data, SDisp_Error);
    //Chi_Square[4] = GetHIVelocity(HIRot_Radii, HIRot_Data, HIRot_Error);
    //Chi_Square[5] = GetHIDispersion(HIDisp_Radii, HIDisp_Data, HIDisp_Error);
    
    //for (int i = 0; i < VC_Radii.size(); ++i)
    //{
    //    cout << VC_Radii_Orig[i] << " " << VC_Radii[i] << endl;
    //}
    
    double chi_square = 0;
    int dof = 0;
    
     for (int i = 0; i < Chi_Square.size(); ++i)
     {
         chi_square += Chi_Square[i];
         dof += DOF[i];
     }
    
    //What follows is based on Binney & Dehnen 1998. Change weights to assign 
    //equal weight to each data set rather than each point.
    double avg_dof = double(dof)/data_sets;
    
    for (int i = 0; i < Chi_Square.size(); ++i)
    {
        chi_square += avg_dof/DOF[i]*Chi_Square[i];
        //dof += DOF[i];
        cout << "DOF " << avg_dof << " " << avg_dof/DOF[i] << endl;
    }
    
    if (do_chisq_file_io)
    {
        ofstream chiout("chisquare.new");
        
        for (int i = 0; i < Chi_Square.size(); ++i)
        {
            chiout << Chi_Square[i] << " " << DOF[i] << endl;
        }
    }
    
    cout << "Chi squares:" << endl;
    
    for (int i = 0; i < Chi_Square.size(); ++i)
    {
        cout << Chi_Square[i] << " " << DOF[i] << endl;
    }
    
    cout << chi_square << endl;
       
    //astro_in.close();
    //error_in.close();
     
    return chi_square;
}

double GetSurfaceBrightness(vector<double>& Radii, vector<double>& Data, 
                            vector<double>& Error)
{
    cout << "Now getting surface brightness profile" << endl;
    
    ofstream sb_out;
    
    if (do_chisq_file_io)
    {
        sb_out.open("surfacebrightness", ios::out);
    }
    
    double chi_square = 0;
    double errorf = Error_Factors[0];
    double integral_bulge = 0;
    double integral_disk[disk], disk_light[disk];
    double kpc_to_arcsec = atan(1.0/Astro.distance)*3600*180/PI;//arcsec per kpc
    
    //We get the projected density in units of mass/kpc^2. This needs to be
    //converted to light using the ML ratio and mass scale (2.325e9 L_solar),
    //to light/arcsec^2 using the above conversion, and then plug into the 
    //equation for magnitude using the distance.
    
    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i];//cout<<errorf<<endl;
        if (radius == 0) continue;
        double error = sqrt(Error[i]*Error[i] + errorf*errorf);
        double point_inclination = acos(1.0-Ellipticity[i]);
                
        //double a_ellip = Radii[i] * sin(Astro.inclination);
        //double b_ellip = sqrt(radius*radius - a_ellip*a_ellip);
        double d_phi = PI/2/99;
        
        for (int j = 0; j < disk; ++j)
        {
            integral_disk[j] = 0;
            disk_light[j] = 0;
        }
    
        for (int j = 0; j < 100; ++j)
        {
            double yp = radius*cos(j*d_phi);
            double xp = radius*sin(j*d_phi)*cos(point_inclination);
            //cout << "step " << i << " " << j << " " << xp << " " << yp << endl;
            
            for (int k = 0; k < disk; ++k)
            {
                double disk_surf_den = DiskSurfaceDensity(xp, yp, k);
                integral_disk[k] += disk_surf_den;
                //cout << integral_disk[k] << endl;
            }
            
            double bulge_surf_den = BulgeSurfaceDensity(xp, yp);
            
            if (j == 0 || j == 99)
            {
                bulge_surf_den *= 0.5;
            }
            
            integral_bulge += bulge_surf_den;
            //cout << "intbulge " << bulge_surf_den << endl;
        }
        
        double bulge_light = 0, disk_light_tot = 0;
        
        for (int j = 0; j < disk; ++j)
        {
            integral_disk[j] /= 100;
            disk_light[j] += integral_disk[j]/Astro.ML_Disk[j]*m_scale/kpc_to_arcsec/kpc_to_arcsec;
            disk_light_tot += disk_light[j];
        }
        
        integral_bulge /= 100;
        bulge_light = integral_bulge/Astro.ml_bulge*m_scale/kpc_to_arcsec/kpc_to_arcsec;
        
        //cout << disk_light_tot << " " << bulge_light << " " << integral_bulge << endl;
        
        //double disk_magnitude_tot = mag_sun+36.572-2.5*log10(disk_light_tot);
        //double bulge_magnitude = mag_sun+36.572-2.5*log10(bulge_light);
        //double total_magnitude = mag_sun+36.572-2.5*log10(disk_light_tot+bulge_light);
        double disk_magnitude_tot = mag_sun-5-2.5*log10(disk_light_tot/Astro.distance/Astro.distance*1e-6);//to 1e-6 is for conv to pc in log arg
        double bulge_magnitude = mag_sun-5-2.5*log10(bulge_light/Astro.distance/Astro.distance*1e-6);
        double total_magnitude = mag_sun-5-2.5*log10((disk_light_tot+bulge_light)/Astro.distance/Astro.distance*1e-6);
        //double disk_magnitude_tot = -27.492-2.5*log10(disk_light_tot/Astro.distance/Astro.distance/206265000/206265000);
        //double bulge_magnitude = -27.492-2.5*log10(bulge_light/Astro.distance/Astro.distance/206265000/206265000);
        //double total_magnitude = -27.492-2.5*log10((disk_light_tot+bulge_light)/Astro.distance/Astro.distance/206265000/206265000);
        double disk_magnitude[disk];
        
        for (int j = 0; j < disk; ++j)
        {
            //disk_magnitude[j] = mag_sun+36.572-2.5*log10(disk_light[j]);
            disk_magnitude[j] = mag_sun-5-2.5*log10(disk_light[j]/Astro.distance/Astro.distance*1e-6);
            //disk_magnitude[j] = -27.492-2.5*log10(disk_light[j]/Astro.distance/Astro.distance/206265000/206265000);
        }

        DOF[0] += 1;
        chi_square += (total_magnitude-Data[i])*(total_magnitude-Data[i])/error/error +
                      2*log(error);
        //cout << "chisq " << chi_square << " " << total_magnitude << " " 
        //     << Radii[i] << " " <<  Data[i] << " " << error << "       " 
        //     << disk_light_tot << " " << bulge_light << endl;
        
        if (do_chisq_file_io)
        {
            sb_out << Radii[i] << " " << Data[i] << " " << Error[i] << " " << error
                   << " " << total_magnitude << " " << disk_magnitude_tot << " ";
            
            for (int j = 0; j < disk; ++j)
            {
                sb_out << disk_magnitude[j] << " ";
            }
                 
            sb_out << bulge_magnitude << endl;
        }
        
        //cout << chi_square << endl;
    }
    
    return chi_square;
}

double GetCircVelocity(vector<double>& Radii, vector<double>& Data, 
                       vector<double>& Error)
{
    cout << "Now getting circular velocities..." << endl;
    
    ofstream vc_out;
    
    if (do_chisq_file_io)
    {
        vc_out.open("circularvelocity", ios::out);
    }
    
    double chi_square = 0;
    double errorf = Error_Factors[1];

    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i];
        double error = sqrt(Error[i]*Error[i] + errorf*errorf);
        //cout << Radii[i]<<" "<<Data[i]<<" "<<Error[i]<<endl;
        double f_rad, f_z, z = 0;
        
        Force(radius, z, f_rad, f_z);
        
        double v_circ = sqrt(-radius*f_rad)*sin_inclination;
        
        DOF[1] += 1;
        chi_square += (v_circ-Data[i])*(v_circ-Data[i])/error/error + 2*log(error);
        //cout << "chisq " << chi_square << " " << v_circ << " " << Radii[i] << " " << Data[i]
        //     << " " << error << endl;
        
        if (do_chisq_file_io)
        {
            vc_out << Radii[i] << " " << v_circ << " " << Data[i] << " " << Error[i] 
                   << " " << error << " " << v_circ/sin_inclination << endl;
        }
        
        //cout << chi_square << endl;
    }
    
    return chi_square;
}

double GetStarVelocity(vector<double>& Radii, vector<double>& Data,
                       vector<double>& Error)
{
    cout << "Now getting stellar line-of-sight velocities..." << endl;
    
    ofstream starvel_out;
    
    if (do_chisq_file_io)
    {
        starvel_out.open("starvelocity", ios::out);
    }
    
    double x = 0, errorf = Error_Factors[2];
    double chi_square = 0;

    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i];
        double data = Data[Data.size()-i-1];
        double error = sqrt(Error[Data.size()-i-1]*Error[Data.size()-i-1] + errorf*errorf);
        //cout << Radii[i]<<" "<<Data[i]<<" "<<Error[i]<<endl;
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, disk_dispersion = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < disk; ++j)
        {
            disk_surface_density1 = DiskSurfaceDensity(x, radius, j);
            disk_surface_density += disk_surface_density1;
            surface_light += disk_surface_density/Astro.ML_Disk[j];
            
            DiskDispersion(x, radius, wbar_disk1, disk_dispersion1, j);
            wbar_disk += wbar_disk1/Astro.ML_Disk[j]*disk_surface_density1;
            //cout << "wbar " << surface_light << " " << wbar_disk << endl;
        } 
        
        //Add on the bulge contribution
        double bulge_surface_density = 0;
        double wbar_bulge = 0, bulge_dispersion;
        
        if (bulge_flag)
        {
            bulge_surface_density = BulgeSurfaceDensity(x, radius);
        
            surface_light += bulge_surface_density/Astro.ml_bulge;
                
            BulgeDispersion(x, radius, wbar_bulge, bulge_dispersion);
            //cout << "wbar " << surface_light << " " << wbar_bulge << endl;
        }
        
        wbar = (wbar_disk + wbar_bulge/Astro.ml_bulge*bulge_surface_density)/surface_light;
        
        DOF[2] += 1;
        //chi_square += (wbar-Data[i])*(wbar-Data[i])/error/error + 2*log(error);
        chi_square += (wbar-data)*(wbar-data)/error/error + 2*log(error);
        
        //cout << "wbar "  << wbar_disk << " " << wbar_bulge << " " << wbar << endl;
        
        //cout << "chisq " << chi_square << " " << wbar << " " << Radii[i] << " " << Data[i]
        //     << " " << data << " " << error << endl;
        
        if (do_chisq_file_io)
        {
            starvel_out << Radii[i] << " " << wbar << " " << wbar_disk1 << " " << wbar_bulge 
                        << " " << data << " " << Error[i] << " " << error << endl;
        }
    }
    
    return chi_square;
}
        
double GetStarDispersion(vector<double>& Radii, vector<double>& Data, 
                         vector<double>& Error)
{
    cout << "Now getting stellar line-of-sight dispersion..." << endl;
    
    ofstream stardisp_out;
    
    if (do_chisq_file_io)
    {
        stardisp_out.open("stardispersion", ios::out);
    }
    
    double chi_square = 0;
    double errorf = Error_Factors[3];

    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i], x = 0;
        double data = Data[Data.size()-i-1];
        double error = sqrt(Error[Data.size()-i-1]*Error[Data.size()-i-1] + errorf*errorf);
        //cout << Radii[i]<<" "<<Data[i]<<" "<<Error[i]<<endl;
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, disk_dispersion = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < disk; ++j)
        {
            disk_surface_density1 = DiskSurfaceDensity(x, radius, j);
            disk_surface_density += disk_surface_density1;
            surface_light += disk_surface_density/Astro.ML_Disk[j];
            
            DiskDispersion(x, radius, wbar_disk1, disk_dispersion1, j);
            wbar_disk += wbar_disk1/Astro.ML_Disk[j]*disk_surface_density1;
            disk_dispersion += disk_dispersion1/Astro.ML_Disk[j]*disk_surface_density1;
        } 
        
        //Add on the bulge contribution
        double bulge_surface_density = 0;
        double wbar_bulge = 0, bulge_dispersion = 0, total_dispersion;
        double disk_surface_light = surface_light, bulge_surface_light = 0;
        
        if (bulge_flag)
        {
            bulge_surface_density = BulgeSurfaceDensity(x, radius);
        
            bulge_surface_light = bulge_surface_density/Astro.ml_bulge;
            surface_light += bulge_surface_light;
                
            BulgeDispersion(x, radius, wbar_bulge, bulge_dispersion);
            //cout << "wbar " << surface_light << " " << wbar_bulge << " " << bulge_dispersion << endl;
        }
        
        wbar = (wbar_disk + wbar_bulge/Astro.ml_bulge*bulge_surface_density)/surface_light;
        
        total_dispersion = (disk_dispersion + bulge_dispersion/Astro.ml_bulge*bulge_surface_density)/
                           surface_light;
                           
        total_dispersion = sqrt(total_dispersion - wbar*wbar);
        
        disk_dispersion = sqrt(disk_dispersion/disk_surface_light - 
                               wbar_disk*wbar_disk/disk_surface_light/disk_surface_light);
        
        if (bulge_flag)
        {
            //cout << "bdisp  " << bulge_dispersion << " " << bulge_surface_light << " "
            //     << wbar_bulge*wbar_bulge/bulge_surface_light/bulge_surface_light << " ";
            
            bulge_dispersion = sqrt(bulge_dispersion - wbar_bulge*wbar_bulge);
                                    //wbar_bulge*wbar_bulge/bulge_surface_light/bulge_surface_light); 
            
            //cout << bulge_dispersion << " " << wbar_bulge << endl;
        }
        
        DOF[3] += 1;
        chi_square += (total_dispersion-data)*(total_dispersion-data)/error/error + 2*log(error);
        
        //cout << "chisq " << chi_square << " " << total_dispersion << " " << disk_dispersion << " " 
        //     << bulge_dispersion << " " << Radii[i] << " " << Data[i] << " " << error << endl;
        
        if (do_chisq_file_io)
        {
            stardisp_out << Radii[i] << " " << total_dispersion << " " << disk_dispersion << " " 
                         << bulge_dispersion << " " << data << " " << Error[i] << " " << error << endl;
        }
        
        //cout << chi_square << endl;
    }
    
    return chi_square;
}

double GetHIVelocity(vector<double>& Radii, vector<double>& Data,
                     vector<double>& Error)
{
    cout << "Now getting HI line-of-sight velocities..." << endl;
    
    ofstream starvel_out;
    
    if (do_chisq_file_io)
    {
        starvel_out.open("HIvelocity", ios::out);
    }
    
    double x = 0, errorf = Error_Factors[4];
    double chi_square = 0;

    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i];
        double data = -Data[i];
        double error = sqrt(Error[i]*Error[i] + errorf*errorf);
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, disk_dispersion = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < disk; ++j)
        {
            disk_surface_density1 = GasDiskSurfaceDensity(x, radius, j);
            disk_surface_density += disk_surface_density1;
            //surface_light += disk_surface_density/Astro.ML_Disk[j];
            
            GasDiskDispersion(x, radius, wbar_disk1, disk_dispersion1, j);
            wbar_disk += wbar_disk1*disk_surface_density1;
        } 
        
        wbar = wbar_disk/disk_surface_density;
        
        DOF[4] += 1;
        //chi_square += (wbar-Data[i])*(wbar-Data[i])/error/error + 2*log(error);
        chi_square += (wbar-data)*(wbar-data)/error/error + 2*log(error);
        
        //cout << "wbar "  << wbar_disk << " " << wbar_bulge << " " << wbar << endl;
        
        //cout << "chisq " << chi_square << " " << wbar << " " << Radii[i] << " " << Data[i]
        //     << " " << data << " " << error << endl;
        
        if (do_chisq_file_io)
        {
            starvel_out << Radii[i] << " " << wbar << " "  
                        << Data[i] << " " << Error[i] << " " << error << endl;
        }
    }
    
    return chi_square;
}
        
double GetHIDispersion(vector<double>& Radii, vector<double>& Data, 
                       vector<double>& Error)
{
    cout << "Now getting HI line-of-sight dispersion..." << endl;
    
    ofstream stardisp_out;
    
    if (do_chisq_file_io)
    {
        stardisp_out.open("HIdispersion", ios::out);
    }
    
    double chi_square = 0;
    double errorf = Error_Factors[5];

    for (int i = 0; i < Data.size(); ++i)
    {
        double radius = Radii[i], x = 0;
        double error = sqrt(Error[i]*Error[i] + errorf*errorf);
        double disk_surface_density = 0, surface_light = 0, disk_surface_density1;
        double wbar_disk = 0, disk_dispersion = 0, wbar = 0, wbar_disk1, disk_dispersion1;
        
        //Add up the contributions to the light and velocities from each disk
        for (int j = 0; j < gas_disk; ++j)
        {
            disk_surface_density1 = GasDiskSurfaceDensity(x, radius, j);
            disk_surface_density += disk_surface_density1;
            //surface_light += disk_surface_density/Astro.ML_Disk[j];
            
            GasDiskDispersion(x, radius, wbar_disk1, disk_dispersion1, j);
            wbar_disk += wbar_disk1*disk_surface_density1;
            disk_dispersion += disk_dispersion1*disk_surface_density1;
        } 
        
        wbar = wbar_disk/disk_surface_density;
        
        double total_dispersion = disk_dispersion/disk_surface_density;
                           
        total_dispersion = sqrt(total_dispersion - wbar*wbar);
        
        DOF[4] += 1;
        chi_square += (total_dispersion-Data[i])*(total_dispersion-Data[i])/error/error + 2*log(error);
        
        //cout << "chisq " << chi_square << " " << total_dispersion << " " << disk_dispersion << " " 
        //     << bulge_dispersion << " " << Radii[i] << " " << Data[i] << " " << error << endl;
        
        if (do_chisq_file_io)
        {
            stardisp_out << Radii[i] << " " << total_dispersion << " " 
                         << disk_dispersion << " " << Data[i] << " " 
                         << Error[i] << " " << error << endl;
        }
        
        //cout << chi_square << endl;
    }
    
    return chi_square;
}

void RotateCoordinates(double &xp, double &yp, double &zp, 
                       double &x, double &y, double &z)
{
    x = xp*cos_inclination - zp*sin_inclination;
    y = yp;
    z = xp*sin_inclination + zp*cos_inclination;
}    
	
void RotateCoordinatesBack(double &xp, double &yp, double &zp, 
                           double &x, double &y, double &z)
{
    x = xp*cos_inclination + zp*sin_inclination;
    y = yp;
    z = -xp*sin_inclination + zp*cos_inclination;
}    

void RotateCoordinatesOffCentre(double &xp, double &yp, double &zp, 
                                double &x, double &y, double &z)
{
    x = xp/cos_inclination - zp*sin_inclination;
    y = yp;
    z = zp*cos_inclination;
}    


// double DiskSurfaceDensity(double &xp, double &yp, int &j)
// {
//     double x, y, z;
//     double sum = 0, zdc = G.Z_Disk[j]/cos_inclination;
//     double t = 2.0/n_los-1, dt = 2*zdc/n_los;
//     double zp = 0.5*zdc*log((1+t)/(1-t));
//     double weights = 0;
//     
//     //RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
//     
//     double r_cyl, disk_dens;
//     //sum += 1.5*disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
//     
//     for (int i = 2; i < n_los-1; ++i)
//     {
//         t = 2.0*i/n_los - 1;
//         zp = 0.5*zdc*log((1+t)/(1-t));
//         //cout << "zp " << i << " " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;
//     
//         RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
//     
//         r_cyl = sqrt(x*x+y*y);
//         disk_dens = DiskDensfI(r_cyl, z, j);
//         sum += disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
//         weights += cosh(zp/zdc)*cosh(zp/zdc);
//         //cout << "zp " << xp << " " << yp << " " << zp << "     " << x << " " 
//         //     << y << " " << z << " " << disk_dens << endl;
//     }
//     
//     //Now do the i=1 and i=n_los-1 cases
//     t = -2.0/n_los + 1;
//     zp = 0.5*zdc*log((1+t)/(1-t));
//     //cout << "zp     " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;
// 
//     RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
// 
//     r_cyl = sqrt(x*x+y*y);
//     disk_dens = DiskDensfI(r_cyl, z, j);
//     sum += 1.5*disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
//     weights += 1.5*cosh(zp/zdc)*cosh(zp/zdc);
//     
//     t = 2.0/n_los - 1;
//     zp = 0.5*zdc*log((1+t)/(1-t));
//     //cout << "zp     " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;
// 
//     RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
//  
//     r_cyl = sqrt(x*x+y*y);
//     disk_dens = DiskDensfI(r_cyl, z, j);
//     sum += 1.5*disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
//     weights += 1.5*cosh(zp/zdc)*cosh(zp/zdc);
//     
//     cout << setprecision(8) << "sum " << sum*dt << " " << xp << " " << yp << endl;
//     
//     return sum*dt;
// }

//This subroutine and the one above give the same answer for the total
//projected mass density in the line of sight, to within about 1-2%. But
//I can't figure out the one above, so I'm using the one below.
double DiskSurfaceDensity(double &xp, double &yp, int &j)
{
    double x, y, z;
    double sum = 0, zdc = G.Z_Disk[j]/cos_inclination;
    double dt = 8*zdc/n_los;
    double zp;
    double weights = 0;
    
    //RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
    
    double r_cyl, disk_dens;
    //sum += 1.5*disk_dens*cosh(zp/zdc)*cosh(zp/zdc);
    
    for (int i = 1; i < n_los; ++i)
    {
        zp = -4*zdc+i*dt;;
        //cout << "zp " << i << " " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;
    
        RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
    
        r_cyl = sqrt(x*x+y*y);
        disk_dens = DiskDensfI(r_cyl, z, j);
        sum += disk_dens*dt;
        //weights += dt;
        //cout << "zp " << xp << " " << yp << " " << zp << "     " << x << " " 
        //     << y << " " << z << " " << disk_dens << " " << dt << endl;
    }
    
    //Now do the i=1 and i=n_los-1 cases
    zp = -4*zdc;
    //cout << "zp     " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;

    RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);

    r_cyl = sqrt(x*x+y*y);
    disk_dens = DiskDensfI(r_cyl, z, j);
    sum += 0.5*disk_dens*dt;
    //weights += 0.5*cosh(zp/zdc)*cosh(zp/zdc);
    
    zp = 4*zdc;
    //cout << "zp     " << zp << " " << zdc << " " << cosh(zp/zdc) << endl;

    RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);

    r_cyl = sqrt(x*x+y*y);
    disk_dens = DiskDensfI(r_cyl, z, j);
    sum += 0.5*disk_dens*dt;
    //weights += 1.5*cosh(zp/zdc)*cosh(zp/zdc);
    
    //cout << setprecision(8) << "sum " << sum << " " << xp << " " << yp << endl;
    
    return sum  ;
}

//This subroutine and the one above give the same answer for the total
//projected mass density in the line of sight, to within about 1-2%. But
//I can't figure out the one above, so I'm using the one below.
double GasDiskSurfaceDensity(double &xp, double &yp, int &j)
{
    double x, y, z;
    double sum = 0, zdc = G.Z_GasDisk[j]/cos_inclination;
    double dt = 8*zdc/n_los;
    double zp;
    double weights = 0;
    
    double r_cyl, disk_dens;
    
    for (int i = 1; i < n_los; ++i)
    {
        zp = -4*zdc+i*dt;;
    
        RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
    
        r_cyl = sqrt(x*x+y*y);
        disk_dens = GasDiskDensfI(r_cyl, z, j);
        sum += disk_dens*dt;
        //weights += dt;
    }
    
    //Now do the i=1 and i=n_los-1 cases
    zp = -4*zdc;

    RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);

    r_cyl = sqrt(x*x+y*y);
    disk_dens = GasDiskDensfI(r_cyl, z, j);
    sum += 0.5*disk_dens*dt;
    //weights += 0.5*cosh(zp/zdc)*cosh(zp/zdc);
    
    zp = 4*zdc;

    RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);

    r_cyl = sqrt(x*x+y*y);
    disk_dens = GasDiskDensfI(r_cyl, z, j);
    sum += 0.5*disk_dens*dt;
    //weights += 1.5*cosh(zp/zdc)*cosh(zp/zdc);
    
    return sum  ;
}

double BulgeSurfaceDensity(double &xp, double &yp)
{
    double x, y, z;
    double sum = 0;
    double rp = xp*xp+yp*yp;
    double rlogmax = log(max(G.a_bulge*pow(10, G.n_sersic), r_max)), rlogmin;
    
    if (rp > 0)
    {
        rlogmin = log(rp/100);
    }
    else 
    {
        rlogmin = log(0.1*dr);
    }
    
    double drlog = (rlogmax-rlogmin)/n_los;
    double deltar = exp(rlogmin);
    
    double zp = 0;
    
    RotateCoordinates(xp, yp, zp, x, y, z);
    
    double r_cyl = sqrt(x*x+y*y);
    double psi = Pot(r_cyl, z);
    double bulge_dens = BulgeDensPsi(psi);
    //cout << "bdens " << psi << " " << bulge_dens << endl;
    sum += 0.5*deltar*bulge_dens;
    
    for (int i = 0; i < n_los; ++i)
    {
        double zlog = rlogmin + i*drlog;
        zp = exp(zlog);
        RotateCoordinatesOffCentre(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        psi = Pot(r_cyl, z);
        bulge_dens = BulgeDensPsi(psi);//cout << "bdens " << psi << " " << bulge_dens << endl;
        sum += drlog*zp*bulge_dens;
    }
    
    return 2*sum;
}

void PreDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, 
                       double &f_max, double &v_circ, int &j)
{
    double freq_omega, freq_kappa, x = pos[0], y = pos[1], z = pos[2];
    double r_cyl = sqrt(x*x+y*y);
    double fs, fz;
    GetOmegaKappa(r_cyl, freq_omega, freq_kappa);
    
    v_phi_max = r_cyl*freq_omega;
    v_cyl[0] = sqrt(SigR2(r_cyl, j));
    v_cyl[1] = 0.5*freq_kappa*v_cyl[0]/freq_omega;
    v_cyl[2] = sqrt(SigZ2(r_cyl, j));
    
    double v_phi_max_old = v_phi_max;
    
    FindMax(r_cyl, z, v_cyl[1], v_phi_max, f_max, j);
    f_max *= 1.1;
    
    v_phi_max = v_phi_max_old;
    
    Force(r_cyl, z, fs, fz);
    v_circ = sqrt(-fs*r_cyl);
    
    //cout << v_phi_max << " " << r_cyl << " " << freq_omega << endl;
}
    
void PreGasDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, 
                          double &f_max, double &v_circ, int &j)
{
    double freq_omega, freq_kappa, x = pos[0], y = pos[1], z = pos[2];
    double r_cyl = sqrt(x*x+y*y);
    double fs, fz;
    GetOmegaKappa(r_cyl, freq_omega, freq_kappa);
    
    v_phi_max = r_cyl*freq_omega;
    v_cyl[0] = sqrt(SigR2Gas(r_cyl, j));
    v_cyl[1] = 0.5*freq_kappa*v_cyl[0]/freq_omega;
    v_cyl[2] = sqrt(SigZ2Gas(r_cyl, j));
    
    double v_phi_max_old = v_phi_max;
    
    FindMaxGas(r_cyl, z, v_cyl[1], v_phi_max, f_max, j);
    f_max *= 1.1;
    
    v_phi_max = v_phi_max_old;
    
    Force(r_cyl, z, fs, fz);
    v_circ = sqrt(-fs*r_cyl);
}
    
void DiskVelocities(double *pos, double &v_phi_max, double *v_cyl, double &f_max,
                    double &v_circ, double &w_los, double &df, int &j)
{
    double g2 = 2, gr, gp, gz;
    
    //cout << "dv " << endl;
    
    while (g2 > 1)
    {
        gr = 8*(gsl_rng_uniform(rand_gen)-0.5);
        gp = 16*(gsl_rng_uniform(rand_gen)-0.5);
        gz = 8*(gsl_rng_uniform(rand_gen)-0.5);
        g2 = gr*gr/16+gp*gp/64+gz*gz/16;//cout << g2 << endl;
    }
        
    //cout << "dv2 " << endl;
    
    double vr = v_cyl[0]*gr;
    double vp = v_phi_max + v_cyl[1]*gp;
    double vz = v_cyl[2]*gz;
    
    //The original code defined FDisk function that simply called DiskDF5ez.
    //So I'm skipping it.
    double r_cyl = sqrt(pos[0]*pos[0]+pos[1]*pos[1]), z = pos[2];
    df = DiskDF5ez(vr, vp, vz, r_cyl, z, j)*v_cyl[0]*v_cyl[1]*v_cyl[2];
    double u, v, u_prime, v_prime;
    
    if (r_cyl != 0)
    {
        u = -vr*pos[0]/r_cyl + vp*pos[1]/r_cyl;
        v = -vr*pos[1]/r_cyl - vp*pos[0]/r_cyl;
    }
    else
    {
        u = (vp-vr)*0.5;
        v = -(vp+vr)*0.5;
    }
    
    RotateCoordinates(u, v, vz, u_prime, v_prime, w_los);
}

void GasDiskVelocities(double *pos, double &v_phi_max, double *v_cyl, double &f_max,
                       double &v_circ, double &w_los, double &df, int &j)
{
    double g2 = 2, gr, gp, gz;
    
    //cout << "dv " << endl;
    
    while (g2 > 1)
    {
        gr = 8*(rand()/float(RAND_MAX)-0.5);
        gp = 16*(rand()/float(RAND_MAX)-0.5);
        gz = 8*(rand()/float(RAND_MAX)-0.5);
        g2 = gr*gr/16+gp*gp/64+gz*gz/16;
    }
        
    double vr = v_cyl[0]*gr;
    double vp = v_phi_max + v_cyl[1]*gp;
    double vz = v_cyl[2]*gz;
    
    //The original code defined FDisk function that simply called DiskDF5ez.
    //So I'm skipping it.
    double r_cyl = sqrt(pos[0]*pos[0]+pos[1]*pos[1]), z = pos[2];
    df = DiskDF5ezGas(vr, vp, vz, r_cyl, z, j)*v_cyl[0]*v_cyl[1]*v_cyl[2];
    double u, v, u_prime, v_prime;
    
    if (r_cyl != 0)
    {
        u = -vr*pos[0]/r_cyl + vp*pos[1]/r_cyl;
        v = -vr*pos[1]/r_cyl - vp*pos[0]/r_cyl;
    }
    else
    {
        u = (vp-vr)*0.5;
        v = -(vp+vr)*0.5;
    }
    
    RotateCoordinates(u, v, vz, u_prime, v_prime, w_los);
}

void DiskDispersion(double &xp, double &yp, double &wbar_disk, 
                    double &disk_dispersion, int &j)
{
    double x, y, z, r_cyl;
    double sum = 0;
    double rp = xp*xp+yp*yp, zdc = G.Z_Disk[j]/cos_inclination;
    double dt = 2*zdc/n_los;
    //double rlogmin, rlogmax = log(G.Z_Disk[j]/cos_inclination);
    
//     if (rp > 0)
//     {
//         rlogmin = log(rp/100);
//     }
//     else 
//     {
//         rlogmin = log(0.1*dr);
//     }
    
    double zoffset = xp*sin_inclination/cos_inclination;
    
    double zp = 0;
    
    wbar_disk = 0;
    disk_dispersion = 0;   
    
    for (int i = 1; i < n_los; ++i)
    {
        double t = 2.0*i/n_los-1;
        zp = zoffset+0.5*zdc*log((1+t)/(1-t));
        
        RotateCoordinates(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);//cout<<xp<<" "<<yp<<" "<<zp<<" "<<x<<" "<<y<<" "<<z<<endl;
        
        double disk_density = DiskDensfI(r_cyl, z, j);
        double fac = cosh(zp/zdc)*cosh(zp/zdc);
        
        double w0 = 0, w1 = 0, w2 = 0, w_los, dfn;
        double v_cyl[3], pos[3] = {x,y,z}, f_max, v_circ, v_phi_max;
        
        PreDiskVelocities(pos, v_phi_max, v_cyl, f_max, v_circ, j);
        
        //cout << "pos " << i << " " << xp << " " << yp << " " << zp << " " 
        //     << x << " " << y << " " << z << endl;
        //cout << "do i " << i << " " << r_cyl << " " << v_phi_max << " " 
        //     << v_cyl[0] << " " << v_cyl[1] << " " << v_cyl[2] << " " << f_max << endl;
        
        for (int k = 0; k < n_vel; ++k)
        {
            DiskVelocities(pos, v_phi_max, v_cyl, f_max, v_circ, w_los, dfn, j);
            w0 += dfn;
            w1 += dfn*w_los;
            w2 += dfn*w_los*w_los;
            //cout << "do " << w_los << " " << dfn << endl;
        }
        
        w1 /= w0;
        w2 /= w0;
        
        sum += disk_density*fac;
        wbar_disk += w1*disk_density*fac;
        disk_dispersion += w2*disk_density*fac;
        
        //cout << "subr " << w0 << " " << w1 << " " << w2 << " " << disk_density 
        //     << " " << fac << " " << cosh(zp/zdc) << " " << zdc << endl;
    }
    
    //cout << "disp " << sum << " " << wbar_disk << " " << disk_dispersion << endl;
    
    wbar_disk /= sum;
    disk_dispersion /= sum;
}

void GasDiskDispersion(double &xp, double &yp, double &wbar_disk, 
                       double &disk_dispersion, int &j)
{
    double x, y, z, r_cyl;
    double sum = 0;
    double rp = xp*xp+yp*yp, zdc = G.Z_GasDisk[j]/cos_inclination;
    double dt = 2*zdc/n_los;
    
    double zoffset = xp*sin_inclination/cos_inclination;
    
    double zp = 0;
    
    wbar_disk = 0;
    disk_dispersion = 0;   
    
    for (int i = 1; i < n_los; ++i)
    {
        double t = 2.0*i/n_los-1;
        zp = zoffset+0.5*zdc*log((1+t)/(1-t));
        
        RotateCoordinates(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        
        double disk_density = GasDiskDensfI(r_cyl, z, j);
        double fac = cosh(zp/zdc)*cosh(zp/zdc);
        
        double w0 = 0, w1 = 0, w2 = 0, w_los, dfn;
        double v_cyl[3], pos[3] = {x,y,z}, f_max, v_circ, v_phi_max;
        
        PreGasDiskVelocities(pos, v_phi_max, v_cyl, f_max, v_circ, j);
        
        for (int k = 0; k < n_vel; ++k)
        {
            GasDiskVelocities(pos, v_phi_max, v_cyl, f_max, v_circ, w_los, dfn, j);
            w0 += dfn;
            w1 += dfn*w_los;
            w2 += dfn*w_los*w_los;
        }
        
        w1 /= w0;
        w2 /= w0;
        
        sum += disk_density*fac;
        wbar_disk += w1*disk_density*fac;
        disk_dispersion += w2*disk_density*fac;
    }
    
    wbar_disk /= sum;
    disk_dispersion /= sum;
}

void BulgeVelocities(double &psi, double &vmag, double &dvmag, double &x, 
                     double &y, double &z, double &w_los, double &dfn)
{
    double vr, vp;
    double r_cyl = sqrt(x*x+y*y);
    double cos_theta = -1+2*(gsl_rng_uniform(rand_gen));
    double r_phi = 2*PI*(gsl_rng_uniform(rand_gen));
    
    double u = vmag*sqrt(1-cos_theta*cos_theta)*cos(r_phi);
    double v = vmag*sqrt(1-cos_theta*cos_theta)*sin(r_phi);
    double w = vmag*cos_theta;//
    
    //cout<<"Rand "<<gsl_rng_uniform(rand_gen)<<endl;
    
    double E = psi - 0.5*vmag*vmag;
    dfn = vmag*vmag*vmag*dvmag*BulgeDF(E);
    //cout << dfn << " " << E << " " << u << " " << v << " " << w << "    " << endl;
    //cout << "vel   " << u*u+v*v+w*w << endl;
    
    if (r_cyl == 0)
    {
        vr = u;
        vp = v;
    }
    else
    {
        vr = (u*x+v*y)/r_cyl;
        vp = (-u*y+v*x)/r_cyl;
    }
    
    if (gsl_rng_uniform(rand_gen) < G.bulge_stream)
    {
        vp = -fabs(vp);
    }
    else
    {
        vp = fabs(vp);
    }
    
    if (r_cyl == 0)
    {
        u = vr;
        v = vp;
    }
    else
    {
        u = (vr*x-vp*y)/r_cyl;
        v = (vr*y+vp*x)/r_cyl;
    }
    
    double uprime, vprime;
    RotateCoordinates(u, v, w, uprime, vprime, w_los);
    //cout << "vel   " << u*u+v*v+w*w << " " << uprime*uprime+vprime*vprime+w_los*w_los << endl;
}

void BulgeDispersion(double &xp, double &yp, double &wbar_bulge, 
                     double &bulge_dispersion)
{
    double x, y, z;
    double sum = 0;
    double rp = sqrt(xp*xp+yp*yp);
    double rlogmax = log(min(G.a_bulge*pow(10, G.n_sersic), r_max)), rlogmin;
  
    if (rp > 0)
    {
        rlogmin = log(rp/100);
    }
    else 
    {
        rlogmin = log(0.1*dr);
    }
    
    double drlog = (rlogmax-rlogmin)/n_los;
    double deltar = exp(rlogmin);
    
    //First integrate over velocity with zp = 0, then zp = rlogmin
    double zp = 0;
    
    RotateCoordinatesBack(xp, yp, zp, x, y, z);
    double r_cyl = sqrt(x*x+y*y);
    double psi = Pot(r_cyl, z);
    double bulge_density = BulgeDensPsi(psi);
    
    double w0 = 0, w1 = 0, w2 = 0, w_los, dfn;
    double v_cyl[3], f_max, v_circ, v_phi_max;
    
    double v_max2 = 2*psi;
    double v_max = sqrt(v_max2), v_min = G.v_bulge/100;
    double vlogmax = log(v_max), vlogmin = log(v_min);
    
    double dvlog = vlogmax - vlogmin;
    
    for (int j = 0; j < n_vel; ++j)
    {
        double vmag = v_min*exp(dvlog*(j+1)/n_vel);
    
        BulgeVelocities(psi, vmag, dvlog, x, y, z, w_los, dfn);
        w0 += dfn;
        w1 += dfn*w_los;
        w2 += dfn*w_los*w_los;
        //cout << j << " " << dfn << " " << w_los*w_los << " " << w0 << " " << w1
        //     << " " << w2 << "   " << psi << " " << vmag << " " << dvlog 
        //     << " " << x << " " << y << " " << z << " " << xp << " " 
        //     << yp << " " << zp << endl;
    }
    
    if (w0 > 0)
    {
        w1 /= w0;
        w2 /= w0;
        sum += deltar*bulge_density;
        wbar_bulge += deltar*bulge_density*w1;
        bulge_dispersion = deltar*bulge_density*w2;
    }
    
    //cout << "       0 " << sum << " " << wbar_bulge << " " 
    //     << bulge_dispersion << " " << deltar*bulge_density << endl;
    
    zp = exp(rlogmin);
    
    RotateCoordinatesBack(xp, yp, zp, x, y, z);
    r_cyl = sqrt(x*x+y*y);
    psi = Pot(r_cyl, z);
    bulge_density = BulgeDensPsi(psi);
    
    w0 = 0;
    w1 = 0;
    w2 = 0;
 
    v_max2 = 2*psi;
    v_max = sqrt(v_max2), v_min = G.v_bulge/100;
    vlogmax = log(v_max), vlogmin = log(v_min);
    
    dvlog = vlogmax - vlogmin;
    
    for (int j = 0; j < n_vel; ++j)
    {
        double vmag = v_min*exp(dvlog*(j+1)/n_vel);
   
        BulgeVelocities(psi, vmag, dvlog, x, y, z, w_los, dfn);
        w0 += dfn;
        w1 += dfn*w_los;
        w2 += dfn*w_los*w_los;
        //cout << j << " " << dfn << " " << w_los*w_los << " " << w0 << " " << w1
        //     << " " << w2 << "   " << psi << " " << vmag << " " << dvlog 
        //     << " " << x << " " << y << " " << z << " " << xp << " " 
        //     << yp << " " << zp <<  endl;
   }
    
    if (w0 > 0)
    {
        w1 /= w0;
        w2 /= w0;
        sum += 0.5*(deltar+drlog*zp)*bulge_density;
        wbar_bulge += 0.5*(deltar+drlog*zp)*bulge_density*w1;
        bulge_dispersion = 0.5*(deltar+drlog*zp)*bulge_density*w2;
    }
    
    //cout << "       1 " << sum << " " << wbar_bulge << " " 
    //     << bulge_dispersion << " " << 0.5*(deltar+drlog*zp)*bulge_density << endl;
    
    //Now do the rest of the points along the LOS
    for (int i = 2; i < n_los; ++i)
    {
        double zlog = rlogmin + i*drlog;
        zp = exp(zlog);
        RotateCoordinatesBack(xp, yp, zp, x, y, z);
        r_cyl = sqrt(x*x+y*y);
        psi = Pot(r_cyl, z);
        bulge_density = BulgeDensPsi(psi);
    
        w0 = 0;
        w1 = 0;
        w2 = 0;
 
        v_max2 = 2*psi;
        v_max = sqrt(v_max2), v_min = G.v_bulge/100;
        vlogmax = log(v_max), vlogmin = log(v_min);
    
        dvlog = vlogmax - vlogmin;
    
        for (int j = 0; j < n_vel; ++j)
        {
            double vmag = v_min*exp(dvlog*(j+1)/n_vel);
        
            BulgeVelocities(psi, vmag, dvlog, x, y, z, w_los, dfn);
            w0 += dfn;
            w1 += dfn*w_los;
            w2 += dfn*w_los*w_los;
            //cout << j << " " << dfn << " " << w_los*w_los << " " << w0 << " " << w1
            //     << " " << w2 << "   " << psi << " " << vmag << " " << dvlog 
            //     << " " << x << " " << y << " " << z << " " << xp << " " 
            //     << yp << " " << zp <<  endl;
        }
    
        //cout << "wbulge " << w0 << " " << w1 << " " << w2 << endl;
        
        if (w0 > 0)
        {
            w1 /= w0;
            w2 /= w0;
            sum += drlog*zp*bulge_density;
            wbar_bulge += drlog*zp*bulge_density*w1;
            bulge_dispersion += drlog*zp*bulge_density*w2;
        }
        
        //cout << "       " << i << " " << sum << " " << wbar_bulge << " " 
        //     << bulge_dispersion << " " << drlog*zp*bulge_density << endl;
    }
    
    wbar_bulge /= sum;
    bulge_dispersion /= sum;
}

void FindMax(double &r, double &z, double &vsigp, double &vpmax, double &fmax, 
             int &j)
{
    double dv = 0.1*vsigp, vpm = vpmax, vpmold, zero=0;
    double v0 = vpm-dv, v1 = vpm+dv;
    double f0 = DiskDF5ez(zero,v0,zero,r,z,j);
    double fmid = DiskDF5ez(zero,vpm,zero,r,z,j);
    double f1 = DiskDF5ez(zero,v1,zero,r,z,j);

    if (fmid>=f0 && fmid>f1)
    {
        fmax = fmid;
    }
    else
    {
        if (f0 > f1)
        {
            double ftmp = f0;
            f0 = f1;
            f1 = ftmp;
            v1 = v0;
            dv = -dv;
        }
        
        vpm = v1;
        
        int flag = 1;
        
        while(flag > 0)
        {
            dv *= 2;
            vpmold = vpm;
            vpm += dv;
            f0 = f1;
            f1 = DiskDF5ez(zero,vpm,zero,r,z,j);
            flag = (f1>f0);
        }
        
        vpmax = vpmold;
        fmax = f0;
    }
}
            
void FindMaxGas(double &r, double &z, double &vsigp, double &vpmax, 
                double &fmax, int &j)
{
    double dv = 0.1*vsigp, vpm = vpmax, vpmold, zero=0;
    double v0 = vpm-dv, v1 = vpm+dv;
    double f0 = DiskDF5ezGas(zero,v0,zero,r,z,j);
    double fmid = DiskDF5ezGas(zero,vpm,zero,r,z,j);
    double f1 = DiskDF5ezGas(zero,v1,zero,r,z,j);

    if (fmid>=f0 && fmid>f1)
    {
        fmax = fmid;
    }
    else
    {
        if (f0 > f1)
        {
            double ftmp = f0;
            f0 = f1;
            f1 = ftmp;
            v1 = v0;
            dv = -dv;
        }
        
        vpm = v1;
        
        int flag = 1;
        
        while(flag > 0)
        {
            dv *= 2;
            vpmold = vpm;
            vpm += dv;
            f0 = f1;
            f1 = DiskDF5ezGas(zero,vpm,zero,r,z,j);
            flag = (f1>f0);
        }
        
        vpmax = vpmold;
        fmax = f0;
    }
}
            
