#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>

using namespace::std;

vector<string> FindParameters(string &line_of_data);

int main(int argc, char *argv[])
{
    int burn_in, burn_in_step, max_correlate_order, paranum;
    char line_of_data[500];
    int string_ceiling = 500;
    double sum_para = 0;
    
    if (argc != 5) 
        {cerr<<"Need the burn-in step, max k, parameter column number and parameter file"<<endl;exit(1);}
    
    burn_in = atoi(argv[1]);
    max_correlate_order = atoi(argv[2]);
    paranum = atoi(argv[3]);
    
    ifstream parameter_file(argv[4], ios::in);
    
    vector<double> para;
    
    while (parameter_file.getline(line_of_data, string_ceiling))
    {
        std::string string_of_data = line_of_data;
        
        //Extract the individual strings delimited by whitespaces
        vector<string> parameter_vector = FindParameters(string_of_data);
        
        burn_in_step = atoi(parameter_vector[0].c_str());
       
        if (burn_in_step < burn_in)
            continue;
        
        //double chi = atof
        
        //Push back vector by the value of parameter
        para.push_back(atof(parameter_vector[paranum].c_str()));
        
        sum_para += atof(parameter_vector[paranum].c_str());
    }
    
    int iterations = para.size();
    
    if (max_correlate_order > iterations) 
        {cerr<<"error: k = "<<max_correlate_order
             <<" larger than iterations = "<<iterations<<endl;exit(1);}
    
    double para_mean = sum_para / iterations;
    
    double variance = 0, covariance = 0, autocorrelation;
    
    cerr << "k-th order autocorrelation:" << endl;
    
    for (int k = 1; k < max_correlate_order; k++)
    {
        for (int t = k; t < iterations; t++)
        {
            covariance += (para[t] - para_mean) * 
                          (para[t-k] - para_mean);
            variance += (para[t] - para_mean) * (para[t] - para_mean);
        }
        
        autocorrelation = covariance / variance;
        
        cout << k << " " << autocorrelation << endl;
    }
        
    
    cerr << "paramean " << para_mean << endl;
    
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
        
    
