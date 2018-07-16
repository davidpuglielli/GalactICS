//Functions to read in parameters and observations. To change the observations
//alter the GetObservedData function.

#include "galaxy.h"
#include "chisq.h"
#include "mcmc.h"

void GetParameters(void)
{
    ifstream dbhin("in.dbh", ios::in);
    ifstream gendenspsiin("in.gendenspsi", ios::in);
    
    std::string string_of_data;
    
    if (!(dbhin.is_open() && gendenspsiin.is_open()))
    {
        cerr << "in.dbh or in.gendenspsi missing. Exiting..." << endl;
        exit(1);
    }
    
    //Read in the halo and bulge parameters from in.dbh
    //The disc parameters are to be read in separately, but the flag remains here
    dbhin >> halo_flag;
    dbhin >> G.c_halo >> G.v_halo >> G.a_halo >> G.drtrunc_halo >> G.cusp 
          >> G.halo_stream;
    dbhin >> disk_flag;
    dbhin >> gasdisk_flag;
    dbhin >> bulge_flag;
    dbhin >> G.n_sersic >> G.ppp >> G.v_bulge >> G.a_bulge >> G.bulge_stream;
    dbhin >> smbh_flag;
    dbhin >> G.bh_mass;
    dbhin >> dr >> nr >> l_max; 
    
    sersic_flag = 0;
    do_file_io = 0;
    do_chisq_file_io = 0;
    nbody_flag = 0;
    chisq_flag = 0;
    contraction_flag = 0;
    contraction_prescription = 0;
    nr_ac = nr;
    
    while (getline(dbhin, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        
        if (input_vector.size()==0) continue;
        
        if (input_vector.at(0)=="sersic_flag")
        {
            sersic_flag = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="do_file_io")
        {
            do_file_io = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="do_chisq_file_io")
        {
            do_chisq_file_io = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="nbody_flag")
        {
            nbody_flag = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="chisq_flag")
        {
            chisq_flag = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="contraction_flag")
        {
            contraction_flag = atoi(input_vector[1].c_str());
        }
        if (input_vector.at(0)=="contraction_prescription")
        {
            contraction_prescription = atoi(input_vector[1].c_str());
            if (contraction_prescription < 0 || contraction_prescription > 3)
            {
	            cout << "Contraction prescription must be 1, 2, or 3.\n";
	            cout << "1 = Blumenthal, 2 = Gnedin, 3 = Abadi" << endl;
                exit(1);
            }
        }
        if (input_vector.at(0)=="nr_ac")
        {
            nr_ac = atoi(input_vector[1].c_str());
            if (nr_ac < nr)
            {
                cout << "nr_ac is less than nr. Resetting to nr..." << endl;
                nr_ac = nr;
            }
        }
    }

    cout << "Read in in.dbh" << endl;

    if (disk_flag == 0 && bulge_flag == 0 && halo_flag == 0 && gasdisk_flag == 0)
    {
        cout << "You need at least one component in in.dbh." << endl;
        exit(1);
    }
    
    //Read in the disk parameters separately. In a different function to
    //allow easy modification with multiple disks.
    if (disk_flag == 1) 
    {
        GetDiskParameters();
    }
        
    //Read in the gas disk parameters separately. In a different function to
    //allow easy modification with multiple disks.
    if (gasdisk_flag == 1) 
    {
        GetGasDiskParameters();
    }
        
    //Get number of points for the various integrations to follow
    //Initialise first
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
    
    while (getline(gendenspsiin, string_of_data))
    {
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        //cerr << input_vector.size() << endl;
        
        if (input_vector.size()==0) continue;
        
        if (input_vector.at(0)=="n_psi")
        {
            n_psi = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_int")
        {
            n_int = atoi(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="max_iter")
        {
            max_iter = atof(input_vector[1].c_str());
        }
        else if (input_vector.at(0)=="n_simpson")
        {
            n_simpson = atof(input_vector[1].c_str());
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
    }
    
    cout << "Read in in.gendenspsi" << endl;
    
    dbhin.close();
    gendenspsiin.close();
}
 
void GetLegendreTable(void)
{
    for (int l = 0; l < 21; l+=2)
    {
        for (int j = 0; j < 1001; ++j)
        {
            double cos_theta = double(j)/1000;
            Legendre_Table[l/2][j] = gsl_sf_legendre_Pl(l, cos_theta);
        }
    }
}
