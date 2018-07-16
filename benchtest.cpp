#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

#define x 100000000

using namespace::std;
double square(int i);

int main(void)
{
    fstream outfile("testbench", ios::out);
    //FILE *out;
    //out = fopen("testbench", "w");
    
    double x2;
    
    for (int i = 0; i < x; ++i)
    {
        //fprintf(out, "%d ", i);
        //outfile << i << " ";
        x2 = square(i);
    }
    
    return 0;
}

double square(int i)
{
    return i*i;
}
