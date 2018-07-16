#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace::std;

vector<string> FindParameters(string &line_of_data);

int main(int argc, char *argv[])
{
    char line_of_data[500];
    int string_ceiling = 500;
    
    int XBin, YBin, Xcolumn, Ycolumn;
    
    if (argc != 6) 
    {
        cerr<<"Need the extracted parameter file xbins, ybins, xcolumn, ycolumn"<<endl;
        exit(1);
    }
    
    //Read in a file produced using extractacceptedparameters and other args
    ifstream parameter_file(argv[1], ios::in);
    XBin = atoi(argv[2]);
    YBin = atoi(argv[3]);
    Xcolumn = atoi(argv[4]);
    Ycolumn = atoi(argv[5]);
    
    //Define a vector for the two parameters we want
    vector<double> Xparameter, Yparameter;
    vector<int> stuck_steps;
    float Contour_Array[XBin][YBin];
    
    double Xmax = 0, Xmin = 1e99, Ymax = 0, Ymin = 1e99;
    
    int count = 0;
    
    while (parameter_file.getline(line_of_data, string_ceiling))
    {
        std::string string_of_data = line_of_data;
        
        //Extract the individual strings delimited by whitespaces
        vector<string> parameter_vector = FindParameters(string_of_data);
        
        //Find maxima for X and Y column data
        if (atof(parameter_vector[Xcolumn-1].c_str()) > Xmax)
            Xmax = atof(parameter_vector[Xcolumn-1].c_str());
        if (atof(parameter_vector[Xcolumn-1].c_str()) < Xmin)
            Xmin = atof(parameter_vector[Xcolumn-1].c_str());
        if (atof(parameter_vector[Ycolumn-1].c_str()) > Ymax)
            Ymax = atof(parameter_vector[Ycolumn-1].c_str());
        if (atof(parameter_vector[Ycolumn-1].c_str()) < Ymin)
            Ymin = atof(parameter_vector[Ycolumn-1].c_str());
        
        Xparameter.push_back(atof(parameter_vector[Xcolumn-1].c_str()));
        Yparameter.push_back(atof(parameter_vector[Ycolumn-1].c_str()));
        stuck_steps.push_back(atoi(parameter_vector[0].c_str()));
        
        count++;
    }
    
    Xmax += 1e-10; Xmin -= 1e-10;
    Ymax += 1e-10; Ymin -= 1e-10;
    
    for (int i = 0; i < XBin; ++i)
    {
        for (int  j= 0; j < YBin; ++j)
        {
            Contour_Array[i][j] = 0;
        }
    }
    
    for (int i = 0; i < count; ++i)
    {
        if (Xparameter[i] >= Xmax || Xparameter[i] <= Xmin ||
            Yparameter[i] >= Ymax || Yparameter[i] <= Ymin)
            continue;
        
        int stuck;
        if (i < count-1) stuck = stuck_steps[i+1] - stuck_steps[i];
        else if (i < count) stuck = 10000 - stuck_steps[i];
        if (stuck==113)cerr<<Xparameter[i]<<" "<<Yparameter[i]<<endl;
        //The X bin 
        int j = (int)((Xparameter[i] - Xmin) / (Xmax - Xmin) * XBin);
        //The Y bin
        int k = (int)((Yparameter[i] - Ymin) / (Ymax - Ymin) * YBin);
        if (stuck==113)cerr<<j<<" "<<k<<endl;
        
        Contour_Array[j][k] += (float)stuck;
    }
    
    for (int i = 0; i < XBin; ++i)
    {
        for (int j = 0; j < YBin; ++j)
        {
            cout << i << " " << j << " " << Contour_Array[i][j] << endl;
        }
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
    
    
