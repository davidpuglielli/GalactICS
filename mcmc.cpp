#include "galaxy.h"
#include "chisq.h"
#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[])
{
	signal(SIGINT, signal_handler);
		
    int total_params, step = 0, accepted_steps = 0, extra_params;
    int continue_flag, chain_length, step_start = 0, pass = 0;
    vector<double> start, sigma, minimum, maximum;
    vector<double> Extra_Parameters; //Q, M_d, M_b, M_h, X
    double old_chisqr, new_chisqr, chitot = 0, chitotsq = 0, scale_factor;
    double temperature;
    //double Chi_Square[data_sets];
    
    std::string string_of_data;
    
    //Get parameter information from file 'mcmc_input'    
    ifstream input_file("mcmc_input", ios::in);
    ifstream input2("mcmc.in", ios::in);
    
    //Get observational data
    //This is the function that needs to be modified if you want more data sets
    //plus the data_sets macro
    //GetObservedData();
        
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
            //cout<<"dpar "<<disk_params<<" "<<start.at(nondisk_params+disk_params)<<endl;
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
            //cout<<"dpar "<<disk_params<<" "<<start.at(nondisk_params+disk_params)<<endl;
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
        else if (mode==6 && input_vector.size()==6)
        {
            //start.push_back(atof(input_vector[2].c_str()));
            //sigma.push_back(atof(input_vector[3].c_str()));
            //minimum.push_back(atof(input_vector[4].c_str()));
            //maximum.push_back(atof(input_vector[5].c_str()));
            
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
    contraction_flag = 0;
    contraction_prescription = 1;
    nr_ac = 10000;

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
        else if (input_vector.at(0)=="contraction_flag")
        {
            contraction_flag = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="contraction_prescription")
        {
            contraction_prescription = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="nr_ac")
        {
            nr_ac = atoi(input_vector[1].c_str());
        }
    }
     
    extra_params = 3+2*disk;
    Extra_Parameters.resize(extra_params);
               
    if (contraction_flag && (contraction_prescription < 0 || contraction_prescription > 3))
    {
	    cout << "Contraction prescription must be 1, 2, or 3.\n";
	    cout << "1 = Blumenthal, 2 = Gnedin, 3 = Abadi" << endl;
        exit(1);
    }

    //Allocate the memory required for all the containers
    AllocateVectors();
    
    //Get observational data
    //This is the function that needs to be modified if you want more data sets
    //plus the data_sets macro
    GetObservedData();
        
    //Get previous parameters if requested; start at last step
    if (continue_flag)
    {
        cout << "Getting old parameter data..." << endl;
    
        //Set position indicator to end of file
        ifstream parafile("parameters.out", ios_base::ate);
        
        if(parafile)
        {
            //Get file size - it's where the position indicator is
            int file_size = parafile.tellg();
            
            //Read backwards until it finds a newline character
            for (int i = file_size-2; i > 0; --i)
            {
                parafile.seekg(i);
                char c = parafile.get();
                
                if (c == '\n')
                    break;
            }
            
            getline(parafile, string_of_data);
            
            vector<string> input_vector = FindParameters(string_of_data);
            
            if (input_vector.size()-4 != start.size()+Extra_Parameters.size())
            {
                cerr << "Number of parameters in input file and parameter file does not\n"
                     << "match. You need the same number of parameters in both." << endl;
                cerr << input_vector.size() << " " << start.size() << " " 
                     << Extra_Parameters.size() << endl;
                exit(1);
            }
            else
            {
                for (int i = 0; i < start.size(); ++i)
                {
                    start.at(i) = atof(input_vector[i+4].c_str());
                }
            }
            
            step_start = atoi(input_vector[0].c_str());
            accepted_steps = atoi(input_vector[1].c_str());
        }
        else
        {
            cerr << "continue_flag is set to true but no parameters.out file found.\n"
                 << "Exiting..." << endl;
            exit(1);
        }
        
        parafile.close();
    }
//     else
//     {
//         ifstream parafile("parameters.out", ios_base::ate);
//         
//         if(parafile)
//         {
//             time_t current_time = time(0);
// 
//             char* date = ctime(&current_time);
// 
//             stringstream new_filename, system_cmnd;
//             new_filename << "parameters.out_" << date << endl;
//             system_cmnd << "mv parameters.out " << 
// 
//             cout << "Old parameters.out renamed to " << new_filename << endl;
//             system("mv pa
    
	//Open the output files, append to parameters.out if needed   
    ofstream parameter_file("parameters.out", ios::app);
    ofstream stat_file("statistics.out", ios::out);
	
    //Get new set of parameters for the first step   
    vector<double> Parameters(total_params), Random(total_params);
    vector<double> Proposal_Parameters(total_params);   
    //double covar[total_params][total_params];
        
    for (int i = 0; i < total_params; ++i)
    {
        Parameters.at(i) = start.at(i);
        Proposal_Parameters.at(i) = Parameters.at(i);
    }
    
    WriteDataFiles(Proposal_Parameters);
        
    double chi_sqr_old = GenGalaxyAndGetChiSquare();
    
    cout << "Got the chisq" << endl;

    //Set up random number generator
    //default is mt19937
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    ////////////////////////////////////////////////////////////////////////
	/////////////////////////// BEGIN MAIN LOOP ////////////////////////////
	////////////////////////////////////////////////////////////////////////
      
    for (step = step_start; step < chain_length; ++step)
    {
        if (step%100 == 0) 
        {
            cout << "Step " << step << " now calculating...  Scale factor is " 
                 << scale_factor << endl;
        }

        //get the set of proposal parameters and check they're in range
        for (int j = 0; j < total_params; ++j)
        {
            Random[j] = gsl_ran_ugaussian(r);
            
            Proposal_Parameters[j] = Parameters[j] + scale_factor*sigma[j]*Random[j];
            
            if (Proposal_Parameters[j] > maximum[j] && sigma[j] != 0)
                Proposal_Parameters[j] = maximum[j];
            else if (Proposal_Parameters[j] < minimum[j] && sigma[j] != 0)
                Proposal_Parameters[j] = minimum[j];
        }
        
        WriteDataFiles(Proposal_Parameters);
        double chi_sqr = GenGalaxyAndGetChiSquare();
        
        //Metropolis-Hastings algorithm
        double diff = chi_sqr - chi_sqr_old;
        double alpha;
        
        if (diff < 0)
        {
            chi_sqr_old = chi_sqr;
            accepted_steps++;
            
            for (int j = 0; j < total_params; ++j)
            {
                Parameters[j] = Proposal_Parameters[j];
                stat_file << Random[j]*Random[j] << " ";
            }
            
            CopyExtraParameters(Extra_Parameters);
        }
        else
        {
            alpha = gsl_rng_uniform(r);
            
            //Use as prior the C-M relation from Maccio et al. - Gaussian fit to distribution of log(c)
//             double mean_concen = GetConcentration(proposal_total_mass);
//             double concentration = proposal_virial_radius/Proposal_Parameters[1];
//             double exparg = (log10(concentration)-log10(mean_concen))*
//                             (log10(concentration)-log10(mean_concen));
//             double prob_concen = exp(-0.5*exparg/sigc/sigc);
//             prob_concen /= sqrt(2*PI)*sigc;
//             
//             double old_mean_concen = GetConcentration(real_masshalo);
//             double old_concentration = real_virial_radius/Parameters[1];
//             double old_exparg = (log10(old_concentration)-log10(old_mean_concen))*
//                                 (log10(old_concentration)-log10(old_mean_concen));
//             double old_prob_concen = exp(-0.5*old_exparg/sigc/sigc);
//             old_prob_concen /= sqrt(2*PI)*sigc;
//              
//             double ratio = prob_concen/old_prob_concen;//cerr << ratio << endl;
            //double ratio = exp(-diff/2/(temperature*exp(step/5000)));
            double ratio = exp(-diff/2);
            double facnew = Proposal_Parameters[25]*Proposal_Parameters[26]*
                            Proposal_Parameters[27]*Proposal_Parameters[28];
            double facold = Parameters[25]*Parameters[26]*
                            Parameters[27]*Parameters[28];
            ratio *= facold/facnew;
            
            if (ratio > alpha)
            {
                chi_sqr_old = chi_sqr;
                accepted_steps++;
                
                for (int i = 0; i < total_params; ++i)
                {
                    Parameters[i] = Proposal_Parameters[i];
                    stat_file << Random[i]*Random[i] << " ";
                }
            
                CopyExtraParameters(Extra_Parameters);
            }
                    
            //stat_file << "rat " << ratio << " " << diff << " " << exp(-diff/2) << " " 
            //          << facold/facnew << " " << alpha << " " << gsl_rng_uniform(rand_gen) << " ";
        }
//                 for (int i = 0; i < total_params; ++i)
//                 {
//                     Parameters[i] = Proposal_Parameters[i];
//                     stat_file << Random[i] << " ";
//                 }

        parameter_file << " " << step << " " << accepted_steps << " " 
                       << chi_sqr_old << " " << chi_sqr << " ";

        for (int i = 0 ; i < total_params; ++i)
        {
            //parameter_file << Proposal_Parameters[i] << " ";
            parameter_file << Parameters[i] << " ";
        }
         
        parameter_file << "    ";
               
        for (int i = 0 ; i < extra_params; ++i)
        {
            //parameter_file << Proposal_Parameters[i] << " ";
            parameter_file << Extra_Parameters.at(i) << " ";
        }
                
        parameter_file << endl;
        
        stat_file << endl;
        
        ClearVectors();
    }
    
    ////////////////////////////////////////////////////////////////////////
	/////////////////////////// END OF MAIN LOOP ///////////////////////////
	////////////////////////////////////////////////////////////////////////
      
    gsl_rng_free(r);
    
    return 0;
}
