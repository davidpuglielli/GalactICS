#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace::std;

vector<string> FindParameters(string &line_of_data);
int compare(const void *a, const void *b);

int main(int argc, char *argv[])
{
    char line_of_data[500];
    int string_ceiling = 500;
    
    int XBin, YBin, Xcolumn, Ycolumn;
    float Xmax, Xmin, Ymax, Ymin;
    
    if (argc != 10) 
    {
        cerr<<"Need the parameter file, xbins, ybins, xcolumn, ycolumn, xmax, xmin, ymax, ymin"<<endl;
        exit(1);
    }
    
    //Read in a file produced using extractacceptedparameters and other args
    ifstream parameter_file(argv[1], ios::in);
    XBin = atoi(argv[2]);
    YBin = atoi(argv[3]);
    Xcolumn = atoi(argv[4]);
    Ycolumn = atoi(argv[5]);
    Xmax = atof(argv[6]);
    Xmin = atof(argv[7]);
    Ymax = atof(argv[8]);
    Ymin = atof(argv[9]);
    
    //Define a vector for the two parameters we want
    vector<double> Xparameter, Yparameter;
    //vector<int> stuck_steps;
    float Contour_Array[XBin][YBin];
    
    //double Xmax = 0, Xmin = 1e99, Ymax = 0, Ymin = 1e99;
    
    int count = 0;
    
    while (parameter_file.getline(line_of_data, string_ceiling))
    {
        std::string string_of_data = line_of_data;
        
        //Extract the individual strings delimited by whitespaces
        vector<string> parameter_vector = FindParameters(string_of_data);
        
        //Find maxima for X and Y column data
        //if (atof(parameter_vector[Xcolumn-1].c_str()) > Xmax)
        //    Xmax = atof(parameter_vector[Xcolumn-1].c_str());
        //if (atof(parameter_vector[Xcolumn-1].c_str()) < Xmin)
        //    Xmin = atof(parameter_vector[Xcolumn-1].c_str());
        //if (atof(parameter_vector[Ycolumn-1].c_str()) > Ymax)
        //    Ymax = atof(parameter_vector[Ycolumn-1].c_str());
        //if (atof(parameter_vector[Ycolumn-1].c_str()) < Ymin)
        //    Ymin = atof(parameter_vector[Ycolumn-1].c_str());
        
        Xparameter.push_back(atof(parameter_vector[Xcolumn-1].c_str()));
        Yparameter.push_back(atof(parameter_vector[Ycolumn-1].c_str()));
        //stuck_steps.push_back(atoi(parameter_vector[0].c_str()));
        
        count++;
    }
    
    cerr << "Xmax " << Xmax << " " << Xmin << endl;
    cerr << "Ymax " << Ymax << " " << Ymin << endl;
    
    //Xmax += 1e-10; Xmin -= 1e-10;
    //Ymax += 1e-10; Ymin -= 1e-10;
    
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
        
        //int stuck;
        //if (i < count-1) stuck = stuck_steps[i+1] - stuck_steps[i];
        //else if (i < count) stuck = 10000 - stuck_steps[i];
        
        //The X bin 
        int j = (int)((Xparameter[i] - Xmin) / (Xmax - Xmin) * XBin);
        //The Y bin
        int k = (int)((Yparameter[i] - Ymin) / (Ymax - Ymin) * YBin);
        
        Contour_Array[j][k] += 1.0;
    }
    
    float unordered[XBin*YBin];
    
    for (int i = 0; i < XBin; ++i)
    {
        for (int j = 0; j < YBin; ++j)
        {
            unordered[i*YBin+j] = Contour_Array[i][j];
            cout << i << " " << j << " " << Contour_Array[i][j] << endl;
//             for (int k = 0; k < i*YBin+j; ++k)
//             {
//                 if (Contour_Array[i][j] < unordered[k])
//                     continue;
//                 else if (Contour_Array[i][j] > unordered[k])
//                 {
//                     int l = i*YBin+j - k;
//                     unordered[l+1] = unordered[l];
//                     unordered[k] = Contour_Array[i][j];
//                     break;
//                 }
//             }
            
            //cout << "          " << i*YBin+j << " " << unordered[i*YBin+j] << endl;
        }
    }
    
    qsort(unordered, XBin*YBin, sizeof(float), compare);
    
    float sum = 0;
    float onesigma = -1, twosigma = -1;
    
    //for (int i = XBin*YBin; i > 0; --i)
    //    cerr<<i<<" "<<unordered[i]<<" "<<count<<endl;
        
    for (int i = XBin*YBin-1; i > -1; --i)
    {
        sum += unordered[i];
        if (sum > 0.6826*count && onesigma==-1)
        {
            onesigma = sum;
            cerr << i << " " << onesigma << " " << unordered[i] << endl;
        }
        
        if (sum > 0.95*count && twosigma==-1)
        {
            twosigma = sum;
            cerr << i << " " << twosigma << " " << unordered[i] << endl;
            break;
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

int compare(const void *a, const void *b)
{
    return ( *(int*)a - *(int*)b );
}
