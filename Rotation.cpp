#include "galaxy.h"

void GetRotationCurve(void)
{
    ofstream rotout("rotation", ios::out);
    
    for (int i = 1; i < Radius.size(); i+=10)
    {
        double f_rad, f_z, z = 0;
        double radius = Radius[i];
        
        Force(radius, z, f_rad, f_z);
        
        double rotation = sqrt(-radius*f_rad);
        
        HaloForce(radius, z, f_rad, f_z);
        
        double halo_rotation = sqrt(-radius*f_rad);
        
        BulgeForce(radius, z, f_rad, f_z);
        
        double bulge_rotation = sqrt(-radius*f_rad);
        
        double disk_rotation = rotation*rotation-halo_rotation*halo_rotation-
                               bulge_rotation*bulge_rotation;
        
        disk_rotation = sqrt(disk_rotation);
        
        rotout << radius << " " << rotation << " " << halo_rotation << " "
               << bulge_rotation << " " < disk_rotation << endl;
    }
}
        
        
        
        
