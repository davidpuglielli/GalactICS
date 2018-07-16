#include <vector>

#define data_sets 4

using namespace::std;

int do_chisq_file_io, n_los, n_vel;

vector<double> Error_Factors(data_sets);
vector<double> Chi_Square(data_sets);
vector<int>    DOF(data_sets);

vector<double> SB_Radii_Orig, P_Angle_Orig, Ellipticity_Orig, SB_Data_Orig, SB_Error_Orig;
vector<double> VC_Radii_Orig, VC_Data_Orig, VC_Error_Orig;
vector<double> SVel_Radii_Orig, SVel_Data_Orig, SVel_Error_Orig;
vector<double> SDisp_Radii_Orig, SDisp_Data_Orig, SDisp_Error_Orig;

vector<double> SB_Radii, P_Angle, Ellipticity, SB_Data, SB_Error;
vector<double> VC_Radii, VC_Data, VC_Error;
vector<double> SVel_Radii, SVel_Data, SVel_Error;
vector<double> SDisp_Radii, SDisp_Data, SDisp_Error;

