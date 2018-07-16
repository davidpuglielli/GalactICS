#include "galaxy.h"
#include "chisq.h"
#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Need model line number and parameterfile, and optionally\n"
             << "mcmc_input file and mcmc.in file" << endl; 
        exit(1);
    }
    
    int total_params, step = 0, accepted_steps = 0, extra_params;
    int continue_flag, chain_length, step_start = 0, pass = 0;
    vector<double> start, sigma, minimum, maximum;
    vector<double> Extra_Parameters; //Q, M_d, M_b, M_h, X
    double old_chisqr, new_chisqr, chitot = 0, chitotsq = 0, scale_factor;
    double temperature;
    
    std::string string_of_data;
    
    int model_no = atoi(argv[1]);
    ifstream parameter_file(argv[2], ios::in);
    ifstream input_file, input2;
    
    //Get parameter information from file 'mcmc_input' if not requested
    //on command line
    if (argc == 4)
    {
        input_file.open(argv[3], ios::in);
        input2.open("mcmc.in", ios::in);
    }
    else if (argc >= 5)
    {
        input_file.open(argv[3], ios::in);
        input2.open(argv[4], ios::in);
    }
    else 
    {
        input_file.open("mcmc_input", ios::in);
        input2.open("mcmc.in", ios::in);
    }
    
    if (!parameter_file || !input_file || !input2)
    {
        cerr << "parameter file, mcmc_input or mcmc.in missing. Exiting... " 
             << endl;
        exit(1);
    } 
    
    int mode = -1;
    nondisk_params = 0;
    disk_params = 0;
    gasdisk_params = 0;
    astro_params = 0;
    error_params = 0;
    
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
            
            mode = 1;//cerr<<"got a halo"<<endl;
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
            
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            disk_params++;
        }
        else if (mode==5 && input_vector.size()==6)
        {
            //Disk_Params is a list of the number of parameters in each disk            
            ++GasDisk_Params[gas_disk-1];
            
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            gasdisk_params++;
        }
        else if (mode==3 && input_vector.size()==6)
        {
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            astro_params++;
        }
        else if (mode==4 && input_vector.size()==6)
        {
            start.push_back(atof(input_vector[2].c_str()));
            sigma.push_back(atof(input_vector[3].c_str()));
            minimum.push_back(atof(input_vector[4].c_str()));
            maximum.push_back(atof(input_vector[5].c_str()));
            
            error_params++;
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
    
    //Get number of points for the various integrations to follow
    //Initialise first
    continue_flag = 0;
    chain_length = 10000;
    temperature = 1;
    scale_factor = 1;
    n_psi = 400;
    n_int = 100;
    max_iter = 20;
    n_simpson = 200;
    fraction = 0.75;
    tolerance_factor = 1e-08;
    n_los = 20;
    n_vel = 200;
    nr_spline = 50;
    n_iter = 10;

    while (getline(input2, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        
        if (input_vector.size()!=2) continue;
        
        if (input_vector.at(0)=="continue_flag")
        {
            continue_flag = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="chain_length")
        {
            chain_length = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="temperature")
        {
            temperature = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="scale_factor")
        {
            scale_factor = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="dr")
        {
            dr = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="nr")
        {
            nr = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="l_max")
        {
            l_max = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="sersic_flag")
        {
            sersic_flag = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="do_file_io")
        {
            do_file_io = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="do_chisq_file_io")
        {
            do_chisq_file_io = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="nbody_flag")
        {
            nbody_flag = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="chisq_flag")
        {
            chisq_flag = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_psi")
        {
            n_psi = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_int")
        {
            n_int = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="max_iter")
        {
            max_iter = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="fraction")
        {
            fraction = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="tolerance_factor")
        {
            tolerance_factor = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_los")
        {
            n_los = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_vel")
        {
            n_vel = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="nr_spline")
        {
            nr_spline = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_iter")
        {
            n_iter = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_simpson")
        {
            n_simpson = atoi(input_vector[1].c_str());
        }
    }
     
    extra_params = 3+2*disk;
    Extra_Parameters.resize(extra_params);
    
    getline(parameter_file, string_of_data);
    
    //Extract the individual strings delimited by whitespaces
    vector<string> input_vector = FindParameters(string_of_data);
    
    if (input_vector.size()-4 != start.size()+Extra_Parameters.size())
    {
        cerr << "Number of parameters in input file and parameter file does not\n"
             << "match. You need the same number of parameters in both." << endl;
        cerr << input_vector.size() << " " << start.size() << " " 
             << Extra_Parameters.size() << endl;
        exit(1);
    }
    
    //get the parameters here
    vector<double> Parameters(start.size());
    
    //Relocate the stream get pointer to the beginning of the file
    //parameter_file.seekg(0, ios::end);
    
    while (getline(parameter_file, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        
        if (atoi(input_vector.at(0).c_str()) < model_no)
        {
            cerr << atoi(input_vector.at(0).c_str()) << endl;
            continue;
        }
        else if (atoi(input_vector.at(0).c_str()) == model_no)
        {
            for (unsigned int j = 0; j < start.size(); ++j)
            {
                Parameters.at(j) = atof(input_vector[j+4].c_str());
            }
            
            break;
        }
        else
        {
            cerr << "No model found. Exiting..." << endl;
            exit(1);
        }
    }
    
    WriteDataFiles(Parameters);
    
    system("./GenerateGalaxy");
    
    return 0;
}
        
    
