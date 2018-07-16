#include "galaxy.h"
#include "chisq.h"
#include "mcmc.h"

vector<string> FindParameters(string &line_of_data)
{
    vector<string> parameterstring;
    typedef string::size_type string_size;
    string_size i = 0;
    
    while (i != line_of_data.size())
    {
        while (i != line_of_data.size() && isspace(line_of_data[i]))
            ++i;
        
        string_size j = i;
        
        while (j != line_of_data.size() && !isspace(line_of_data[j]))
            ++j;
        
        if (i != j)
        {
            parameterstring.push_back(line_of_data.substr(i, j - i));
            i = j;
        }
    }
    
    return parameterstring;
}

double GenGalaxyAndGetChiSquare() 
{
    //Uncomment the srand() call if you want to reset the RNG. This is useful
    //if you want the chi square calculations to be repeatable, since there
    //are multiple calls to rand() in the velocity routines.
    //srand(1);
    
    GetParameters();
    DBH();
    GetFreqs();
    
    if (disk_flag || gasdisk_flag)
    {
        DiskDF();
    }

    double chisqr = GetChiSquare(Chi_Square);
    
    //ClearVectors();

    return chisqr;
}

void CopyExtraParameters(vector<double> &Extra_Params)
{
    Extra_Params.at(0) = halo_mass;
    Extra_Params.at(1) = total_disk_mass;
    Extra_Params.at(2) = bulge_mass;
    
    for (int i = 0; i < disk; ++i)
    {
        Extra_Params.at(3+2*i) = Q_avg.at(i);
        Extra_Params.at(4+2*i) = X_avg.at(i);
    }
}

double GetConcentration(const double mass)
{
    double mass1 = mass*2.325e9;
	double conc = pow(10, zero)*pow(mass1/1e12*hubble, slope);

	return conc;
}

void WriteDataFiles(vector<double> &param)
{
    ofstream data_file("in.dbh", ios::out);
    ofstream data_file2("in.diskpars", ios::out);
    ofstream data_file3("in.astronomical", ios::out);
    ofstream data_file4("in.errorfactors", ios::out);
    ofstream data_file5("in.gendenspsi", ios::out);
    ofstream data_file6("in.gasdiskpars", ios::out);
    
    //output in.dbh
    data_file << halo_flag << endl;
    data_file << param.at(0) << " " << param.at(1) << " " << param.at(2) << " "
              << param.at(3) << " " << param.at(4) << " " << param.at(5) << endl;
    data_file << disk_flag << endl;
    data_file << gasdisk_flag << endl;
    data_file << bulge_flag << endl;
    data_file << param.at(6) << " " << param.at(7) << " " << param.at(8) << " "
              << param.at(9) << " " << param.at(10) << endl;
    data_file << smbh_flag << endl;
    data_file << param.at(11) << endl;
    data_file << dr << " " << nr << " " << l_max << " " << endl;
    data_file << "sersic_flag " << sersic_flag << "\n" 
              << "do_file_io " << do_file_io << "\n" 
              << "do_chisq_file_io " << do_chisq_file_io << "\n" 
              << "nbody_flag " << nbody_flag << "\n" 
              << "chisq_flag " << chisq_flag << "\n"
              << "contraction_flag " << contraction_flag << "\n"
              << "contraction_prescription " << contraction_prescription << "\n"
              << "nr_ac " << nr_ac << endl;
    
    cout << "Wrote in.dbh " << endl;
    
    //output in.diskpars
    int k = nondisk_params;
    
//     for (int i = 0; i < Disk_Params.size(); ++i)
//     {
//         if (i > 0)
//             k += Disk_Params[i-1];
//     
//         for (int j = 0; j < Disk_Params[i]; ++j)
//         {
//             int k1 = k+j;
//             data_file2 << param.at(k1) << " ";
//             cout << j << " " << k1 << " " << param.at(k1) << endl;
//         }
//         
//         data_file2 << 0 << endl;
//     }
   
    //Ensure that R_disk = R_sigma
    double rsigma;
    for (int i = 0; i < Disk_Params.size(); ++i)
    {
        if (i > 0)
            k += Disk_Params[i-1];
    
        for (int j = 0; j < Disk_Params[i]; ++j)
        {
            int k1 = k+j;
            
//             if (j == 1)
//             {
//                 rsigma = param.at(k1);
//             }
            
            data_file2 << param.at(k1) << " ";
            //cout << j << " " << k1 << " " << param.at(k1) << endl;
        }
        
        data_file2 << 0 << endl;
    }
    
    cout << "Wrote in.diskpars" << endl;
    
    k = nondisk_params + disk_params;
    
    //output in.gasdiskpars
    for (int i = 0; i < GasDisk_Params.size(); ++i)
    {
        if (i > 0)
            k += GasDisk_Params[i-1];
    
        for (int j = 0; j < GasDisk_Params[i]; ++j)
        {
            int k1 = k+j;
            data_file6 << param.at(k1) << " ";
            cout << j << " " << k1 << " " << param.at(k1) << endl;
        }
        
        data_file6 << 0 << endl;
    }
   
    cout << "Wrote in.gasdiskpars" << endl;
    
    k = nondisk_params + disk_params + gasdisk_params;
    
    //output in.astronomical
    for (int i = 0; i < astro_params; ++i)
    {
        int k1 = k+i;
        data_file3 << param.at(k1) << endl;
        //cout << i << " " << k1 << endl;
    }
    
    cout << "Wrote in.astronomical" << endl;
    
    k = nondisk_params + disk_params + gasdisk_params + astro_params;
    
    //output in.errors
    for (int i = 0; i < error_params; ++i)
    {
        int k1 = k+i;
        data_file4 << param.at(k1) << endl;
        //cout << i << " " << k1 << endl;
    }
    
    cout << "Wrote in.errorfactors" << endl;
    
    data_file5 << "n_psi " << n_psi << endl
               << "n_int " << n_int << endl
               << "max_iter " << max_iter << endl
               << "n_simpson " << n_simpson << endl
               << "tolerance_factor " << tolerance_factor << endl
               << "n_los " << n_los << endl
               << "n_vel " << n_vel << endl
               << "nr_spline " << nr_spline << endl
               << "n_iter " << n_iter << endl
               << "fraction " << fraction << endl;
    
    data_file.close();
    data_file2.close();
    data_file3.close();
    data_file4.close();
    data_file5.close();
    data_file6.close();
}

void signal_handler(int signum)
{
	cerr << "Caught signal " << signum << endl;
	cerr << "Terminating program" << endl;
	exit(signum);
}
