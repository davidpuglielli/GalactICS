//GetObservedData reads in in.astronomical and the data, converting units if
//necessary. This is the function to modify if you want to change the input
//observations, in addition to the convert calls in GetChiSquare().

#include "galaxy.h"
#include "chisq.h"

void GetObservedData(void)
{
    ifstream gasrotation("NGC2841-HI-THINGS.txt", ios::in);
    ifstream starrotation("NGC2841-starrot", ios::in);
    ifstream surfacebrightness("NGC2841-sb-I", ios::in);
    ifstream stardispersion("NGC2841-stardisp", ios::in);
    //ifstream HIrotation("NGC2841-HI-rotation", ios::in);
    //ifstream HIdispersion("NGC2841-HI-disp", ios::in);
	
    ReadData(gasrotation, VC_Radii_Orig, VC_Data_Orig, VC_Error_Orig);
    ReadData(starrotation, SVel_Radii_Orig, SVel_Data_Orig, SVel_Error_Orig);
    ReadData(surfacebrightness, SB_Radii_Orig, P_Angle, Ellipticity, SB_Data, SB_Error);
    ReadData(stardispersion, SDisp_Radii_Orig, SDisp_Data_Orig, SDisp_Error_Orig);
    //ReadData(HIrotation, HIRot_Radii_Orig, HIRot_Data_Orig, HIRot_Error_Orig);
    //ReadData(HIdispersion, HIDisp_Radii_Orig, HIDisp_Data_Orig, HIDisp_Error_Orig);
}

//Overload the ReadData function name to handle 3-5 columns of data
void ReadData(ifstream& file, vector<double>& Vec_1, 
              vector<double>& Vec_2, vector<double>& Vec_3)
{
    double data;
    int count = 0;
    cout << "here" << endl;
    while (file >> data)
    {
        if (count % 3 == 0) Vec_1.push_back(data);//cout << "here " << count << " " << data << endl;
        if (count % 3 == 1) Vec_2.push_back(data);
        if (count % 3 == 2) Vec_3.push_back(data);
        count++;
    }
}   
	
void ReadData(ifstream& file, vector<double>& Vec_1, vector<double>& Vec_2, 
              vector<double>& Vec_3, vector<double>& Vec_4)
{
    double data;
    int count = 0;
    
    while (file >> data)
    {
        if (count % 4 == 0) Vec_1.push_back(data);
        if (count % 4 == 1) Vec_2.push_back(data);
        if (count % 4 == 2) Vec_3.push_back(data);
        if (count % 4 == 3) Vec_4.push_back(data);
        count++;
    }
}   
	
void ReadData(ifstream& file, vector<double>& Vec_1, vector<double>& Vec_2, 
              vector<double>& Vec_3, vector<double>& Vec_4, vector<double>& Vec_5)
{
    double data;
    int count = 0;
    
    while (file >> data)
    {
        if (count % 5 == 0) Vec_1.push_back(data);
        if (count % 5 == 1) Vec_2.push_back(data);
        if (count % 5 == 2) Vec_3.push_back(data);
        if (count % 5 == 3) Vec_4.push_back(data);
        if (count % 5 == 4) Vec_5.push_back(data);
        count++;
    }
}   
	
double ArcsecToKpc(double &rad)
{
    double radius = Astro.distance*tan(4.848136811e-6*rad);
    return radius;
}

//Assumes r_sys is in arcsec. Otherwise, how do you know what it is in kpc?        
void ConvertRadius(vector<double> &Vec_Orig, vector<double> &Vec)
{
    if (Vec.size() != Vec_Orig.size())
    {
        Vec.resize(Vec_Orig.size());
    }
    
    double rsys = ArcsecToKpc(Astro.r_sys);
    
    for (int i = 0; i < Vec_Orig.size(); ++i)
    {
        Vec[i] = ArcsecToKpc(Vec_Orig[i]);
        Vec[i] -= rsys;
    }
}

//Assumes v_sys is in km/s. Otherwise, be careful to comment out the vsys 
//declaration line below
void ConvertVelocities(vector<double> &Vec1_Orig, vector<double> &Vec2_Orig, 
                       vector<double> &Vec1, vector<double> &Vec2)
{
    if (Vec1.size() != Vec2.size())
    {
        cerr << "Warning: Velocity conversion vectors do not match. "
             << "Make sure you have velocity and\nerror columns of the "
             << "same size in your data files. Continuing..." << endl;
    }
    
    if (Vec1.size() != Vec1_Orig.size() || Vec2.size() != Vec2_Orig.size())
    {
        Vec1.resize(Vec1_Orig.size());
        Vec2.resize(Vec2_Orig.size());
    }
    
    double vsys = Astro.v_sys/100.0;
    
    for (int i = 0; i < Vec1_Orig.size(); ++i)
    {
        Vec1[i] = Vec1_Orig[i]/100.0;
        Vec1[i] -= vsys;
    }
    
    for (int i = 0; i < Vec2.size(); ++i)
    {
        Vec2[i] = Vec2_Orig[i]/100.0;
    }
}
