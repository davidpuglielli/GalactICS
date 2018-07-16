#include "galaxy.h"
#include "chisq.h"
#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 4) 
    {
        cerr << "Need the burn-in step, number of bins, parameter file "
             << "and optionally\nMCMC input file. Exiting..." << endl;
        exit(1);
    }
    
    int burn_in = atoi(argv[1]);
    int n_bin = atoi(argv[2]);
    ifstream parameter_file(argv[3], ios::in);
    ifstream input_file;
    
    if (argc == 4)
    {
        input_file.open("mcmc_input", ios::in);
    }
    else if (argc > 4)
    {
        input_file.open(argv[4], ios::in);
    }
    
    if (!parameter_file || !input_file)
    {
        cerr << "Parameter file or input file missing. Exiting..." << endl;
        exit(1);
    }
    
    int total_params, step = 0, accepted_steps = 0, extra_params;
    int continue_flag, chain_length, step_start = 0, pass = 0;
    vector<double> start, sigma, minimum, maximum;
    //vector<double> Extra_Parameters; //Q, M_d, M_b, M_h, X
    vector<string> parameter_names;
    double old_chisqr, new_chisqr, chitot = 0, chitotsq = 0, scale_factor;
    double temperature;
    
    std::string string_of_data;
    
    int mode = -1;
    nondisk_params = 0;
    disk_params = 0;
    gasdisk_params = 0;
    astro_params = 0;
    error_params = 0;
    extra_params = 0;
    
    //cerr << "fhere" << endl;
    
    while (getline(input_file, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        
        if (input_vector.size()==0) continue;
        
        //Only include the information on the galaxy parameters
        //This big if block finds out if the respective components are included
        //and defines a mode flag to be used in the next block to determine
        //whether to skip or include those parameters, and which components
        //those parameters belong to
		if (input_vector[0]=="Halo")
		{
		    halo_flag = atoi(input_vector[1].c_str());
            
		    //if (halo_flag==0)
            //{
            //    mode = 0;
            //    continue;
            //}
            
            mode = 1;
            continue;
        }
		else if (input_vector[0]=="Bulge")
		{
		    bulge_flag = atoi(input_vector[1].c_str());
            
		    //if (bulge_flag==0)
            //{
            //    mode = 0;
            //    continue;
            //}
            
            mode = 1;
		    continue;
		}
		else if (input_vector[0]=="BlackHole")
		{
		    smbh_flag = atoi(input_vector[1].c_str());
            
		    //if (smbh_flag==0)
            //{
            //    mode = 0;
            //    continue;
            //}
            
            mode = 1;
		    continue;
		}
		else if (input_vector[0]=="Disk")
		{
		    if (atoi(input_vector[1].c_str())==0)
            {
                mode = 0;
                continue;
            }
            
            disk_flag = 1;
            Disk_Params.push_back(0);
            mode = 2;
            ++disk;
		    continue;
		}
		else if (input_vector[0]=="GasDisk")
		{
		    if (atoi(input_vector[1].c_str())==0)
            {
                mode = 0;
                continue;
            }
            
            gasdisk_flag = 1;
            GasDisk_Params.push_back(0);
            mode = 5;
            ++gas_disk;
		    continue;
		}
		else if (input_vector[0]=="Astro")
		{
            mode = 3;
		    continue;
		}
		else if (input_vector[0]=="Error")
		{
            mode = 4;
		    continue;
		}
		else if (input_vector[0]=="Extra")
		{
            mode = 6;
		    continue;
		}
        
        //Now either read the parameters if the flags are true (ie. mode > 0)
        //or skip if the flags are false (ie. mode == 0)
        if (mode==-1)
        {
            cerr << "Input file must start with a component name. See README for info" << endl;
            exit(1);
        }
        else if (mode==0)
        {
            continue;
        }
        else if (mode==1 && input_vector.size()==6)
        {
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            nondisk_params++;
        }
        else if (mode==2 && input_vector.size()==6)
        {
            //Disk_Params is a list of the number of parameters in each disk            
            ++Disk_Params[disk-1];
            
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            disk_params++;
        }
        else if (mode==5 && input_vector.size()==6)
        {
            ++GasDisk_Params[gas_disk-1];
            
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            gasdisk_params++;
        }
        else if (mode==3 && input_vector.size()==6)
        {
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            astro_params++;
        }
        else if (mode==4 && input_vector.size()==6)
        {
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            error_params++;
        }
        else if (mode==6 && input_vector.size()==6)
        {
            parameter_names.push_back(input_vector[0]);
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            extra_params++;
        }
    }

    if (astro_params-5 != Disk_Params.size())
    {
        cerr << "Incorrect number of astro parameters. You need as many\n"
             << "disk M/L parameters in mcmc_input as disks" << endl;
        exit(1);
    }

    if (error_params != data_sets)
    {
        cerr << "Incorrect number of error parameters. You need as many error\n"
             << "parameters in mcmc_input as data sets" << endl;
        exit(1);
    }

    total_params = nondisk_params + disk_params + astro_params + error_params;
    
    //extra_params = 3+2*disk;
    
    if (extra_params != 3+2*disk)
    {
        cerr << "Incorrect number of extra parameters. You need 3+2*(number_of_disks)\n"
             << "extra parameters in mcmc_input" << endl;
        exit(1);
    }

    //Extra_Parameters.resize(extra_params);
    
    //for (int j = 0; j <extra_params
    
    double Parameters[total_params+extra_params][n_bin];
    
    for (int i = 0; i < total_params+extra_params; i++)
        for (int j = 0; j < n_bin; j++)
            Parameters[i][j] = 0;
    
    double Mean[total_params+extra_params], Square[total_params+extra_params];
    double Dispersion[total_params+extra_params];
    
    for (int i = 0; i < total_params+extra_params; i++)
    {
        Mean[i] = Square[i] = 0;
    }
    
    int total_iterations = 0;
    
    while (getline(parameter_file, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> parameter_vector = FindParameters(string_of_data);
        
        int burn_in_step = atoi(parameter_vector[0].c_str());
       
        if (burn_in_step < burn_in)
            continue;
        
        for (int j = 0; j < parameter_vector.size()-4; ++j)
        {
            double parameter = atof(parameter_vector[j+4].c_str());
            
            Mean[j] += parameter;
            Square[j] += parameter*parameter;
            
            int k = (int)((parameter-minimum.at(j))/
                    (maximum.at(j)-minimum.at(j))*n_bin);
            
            if (k > n_bin) 
            {
                cerr << "Out of range! Too high! Parameter " << parameter_names[j]
                     << " " << total_iterations << " " << k << " " << parameter 
                     << "        "  << minimum[j] << " " << maximum[j] << endl; 
                continue;
            }
            if (k < 0) 
            {
                cerr << "Out of range! Too low! Parameter " << parameter_names[j] 
                     << " " << total_iterations << " " << k << " " << parameter
                     << "        "  << minimum[j] << " " << maximum[j] << endl; 
                continue;
            }
            else
//             {
//                 cerr << "In range! Just right! Parameter " << parameter_names[j] 
//                      << " " << total_iterations << " " << k << " " << parameter 
//                      << "        "  << minimum[j] << " " << maximum[j] << endl; 
//             }
            
            Parameters[j][k] += 1.0;
        }
        
        ++total_iterations;
    }
    
    for (int j = 0; j < n_bin; j++)
    {
        cout << j << "   ";
        
        for (int i = 0; i < total_params+extra_params; i++)
        {
            cout << (j+0.5)*(maximum[i]-minimum[i])/n_bin+minimum[i] << " " 
                 << Parameters[i][j]/total_iterations*n_bin/(maximum[i]-minimum[i]) 
                 << "    ";//for prob density of 1
        }
        
        cout << endl;
    }
    
    for (int j = 0; j < total_params+extra_params; ++j)
    {
        Mean[j] /= total_iterations;
        Square[j] /= total_iterations;
        
        double variance = Square[j]-Mean[j]*Mean[j];
        Dispersion[j] = sqrt(variance);
        
        cerr << "$ " << setw(12) << left << parameter_names[j] << " $ & $ " 
             << setw(12) << right << Mean[j] << " $ & $ " << setw(12) 
             << variance << " \\pm " << setw(12) << Dispersion[j] << " $" << endl;
    }
    
    return 0;
}
