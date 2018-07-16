// code snippet to get the frequencies from either the potential/density vectors or
// from dbh.dat

#include "galaxy.h"

void GetFreqs(void)
{
    cout<<"Now running GetFreqs..."<<endl;
    
    if (do_file_io)
    {
        ReadFreqs();
    }
    else
    {
        for (int i = 0; i < nr; ++i)
        {
            float rad = Radius[i];
            
            //Total quantities            
            Pot_Major_Tot[i] = PotVector(rad, 0, A_Pot);//cout<<i<<" "<<Pot_Major_Tot[i]<<" ";
            Pot_Minor[i] = PotVector(0, rad, A_Pot);//cout<<Pot_Minor[i]<<" ";
            VC_Tot[i] = sqrt(max(0.0, rad*PotVector(rad, 0, F_R)));//cout<<VC_Tot[i]<<" ";
            Psi2[i] = PotVector(rad, 0, FR2);//cout<<Psi2[i]<<endl;
            
            //Halo only quantities
            Pot_Major[i] = PotVector(rad, 0, Halo_Pot);//cout<<"   "<<Pot_Major[i]<<" ";
            Pot_Up[i] = PotVector(rad*cos(0.05), rad*sin(0.05), Halo_Pot);//cout<<Pot_Up[i]<<" ";
            VC_Major[i] = sqrt(max(0.0, rad*PotVector(rad, 0, Halo_FR)));//cout<<VC_Major[i]<<" ";
            
            if (i==0)
            {
                Vert_Freq[i] = 0;//cout<<Vert_Freq[i]<<endl;
            }
            else
            {
                Vert_Freq[i] = sqrt(max(0.0, VC_Major[i]*VC_Major[i]+
                               800*(Pot_Up[i]-Pot_Major[i])))/rad;//cout<<Vert_Freq[i]<<endl;
            }
                           
            //Bulge only quantities
            Pot_Major[i] = PotVector(rad, 0, Bulge_Pot);
            Pot_Up[i] = PotVector(rad*cos(0.05), rad*sin(0.05), Bulge_Pot);
            VC_Bulge[i] = sqrt(max(0.0, rad*PotVector(rad, 0, Bulge_FR)));
            
            if (i==0)
            {
                Vert_Freq_Bulge[i] = 0;
            }
            else
            {
                Vert_Freq_Bulge[i] = sqrt(max(0.0, VC_Bulge[i]*VC_Bulge[i]+
                                     800*(Pot_Up[i]-Pot_Major[i])))/rad;
            }
        }
    }
}
        
            
    
    
       
        
        
        
