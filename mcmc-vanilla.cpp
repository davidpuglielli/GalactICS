#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <cctype>
#include <stdlib.h>

using namespace::std;

vector<string> FindParameters(string &line_of_data);
void WriteDataFiles(vector<double> *param);

int main(int argc, char *argv[])
{
    int count = 0, step = 0;
    vector<double> start, sigma, minimum, maximum;
    double random_number = 0.5;    
    
    char line_of_data[500];
    int string_ceiling = 500;
    
    ifstream input_file("mcmc_input", ios::in);
    ofstream parameter_file("parameters.out", ios::out);
    
    //Get parameter information from file 'start'    
    while (input_file.getline(line_of_data, string_ceiling))
    {
        //Copy from char* to string
        std::string string_of_data = line_of_data;
        
        //Extract the individual strings delimited by whitespaces
        vector<string> input_vector = FindParameters(string_of_data);
        
        //Only include the information on the galaxy parameters
        if (input_vector.size() == 1) 
            break;
        
	    start.push_back(atof(input_vector[2].c_str()));
	    sigma.push_back(atof(input_vector[3].c_str()));
	    minimum.push_back(atof(input_vector[4].c_str()));
	    maximum.push_back(atof(input_vector[5].c_str()));
	
        count++;
    }

    input_file >> chainlength;
        
    //write the starting data
    WriteDataFiles(&start); 
 
    double parameters[count], proposal_parameters[count], random[count];
    double covar[count][count];
    
    for (int i = 0; i < count; ++i)
    {
        parameters[i] = start[i];
        proposal_parameters[i] = parameters[i];
    }

    //now do the last pass at temperature one         
    for (int step = 0; step < chainlength; ++step)
    {
        //get a sequence of random numbers
        //GetRandomNumbers();

        //get the set of proposal parameters and check they're in range
        for (int j = 0; j < count; ++j)
        {
            proposal_parameters[j] = proposal_parameters[j] + sigma[j] * (rand()-RAND_MAX/2)/RAND_MAX*2;
            if (proposal_parameters[j] > maximum[j])
                proposal_parameters[j] = maximum[j];
            else if (proposal_parameters[j] < minimum[j])
                proposal_parameters[j] = minimum[j];
                
            parameters[j] = proposal_parameters[j];
        }

        WriteDataFiles(parameters);
        system("./einastoprofilelog");
        system("./einastofit1");

        parameter_file << " " << step << " ";

        for (int i = 0 ; i < count; ++i)
        {
            parameter_file << parameters[i] << " ";
        }

        double r0, rho0, alpha;
        ifstream fitted("fittedpeak", ios::in);

        fitted >> r0 >> rho0 >> alpha;
        
        parameter_file << r0 << " " << rho0 << " " << alpha << endl;
    }
    
    return 0;
}

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

void WriteDataFiles(vector<double> *param)
{
    ofstream data_file("einastoparam.in", ios::out);
    ofstream data_file2("fit.einastoparam", ios::out);

    data_file << param->at(0) << " " << param->at(1) << " " << param->at(2) 
              << " 0.0001 " << param->at(3) << " " << param->at(4) << " " << param->at(5) << endl;

    data_file2 << param->at(0) << " " << param->at(1) << " " << param->at(2)<< endl; 
    
    data_file.close();
    data_file2.close();
}
