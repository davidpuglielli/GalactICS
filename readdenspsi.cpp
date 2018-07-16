#include "galaxy.h"

//Why do we need these files?

void ReadDensPsiBulge(double *TableE, double *DensPsiBulge, double *DFSersic)
{
    FILE *denspsi, *df;
    denspsi = fopen("denspsibulge.dat", "r");
    df = fopen("dfsersic.dat", "r");
    
    for (int i = 0; i < n_psi; ++i)
    {
        fscanf(denspsi, "%lf %lf", &TableE[i], &DensPsiBulge[i]);
        fscanf(df, "%lf %lf", &TableE[i], &DFSersic[i]);
    }
}
        
void ReadDensPsiHalo(double *TableE, double *DensPsiHalo, double *DF_NFW)
{
    FILE *denspsi, *df;
    denspsi = fopen("denspsihalo.dat", "r");
    df = fopen("dfhalo.dat", "r");
    
    for (int i = 0; i < n_psi; ++i)
    {
        fscanf(denspsi, "%lf %lf", &TableE[i], &DensPsiHalo[i]);
        fscanf(df, "%lf %lf", &TableE[i], &DF_NFW[i]);
    }
}
 
