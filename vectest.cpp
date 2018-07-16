#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

#define NN 10000000

using namespace::std;

int main(void)
{
    int count = 0;
    vector<double> apot(NN);
    //apot.reserve(NN);
    cout<<apot.max_size();
    //double *apot = new double[NN];

    //for (int i = 0; i < NN; ++i)
    for (vector<double>::iterator it = apot.begin(); it != apot.end(); ++it)
    {
        //apot.push_back(0.28209479177387814347*sqrt(2*i+1));
        *it = (0.28209479177387814347*sqrt(2*count+1));
        //apot[i] = 0.28209479177387814347*sqrt(2*i+1);
        ++count;
    }

    return 0;
}
