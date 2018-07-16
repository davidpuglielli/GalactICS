#include "galaxy.h"
#include "chisq.h"
#include "nbody.h"

using namespace std;

int main(void)
{
    GetParameters();
    AllocateVectors();
    
    DBH();
    GetFreqs();
    
    if (disk_flag || gasdisk_flag)
    {
        DiskDF();
    }
    
    if (nbody_flag)
    {
        GetNBodyRealisations();
    }
    
    if (chisq_flag)
    {
        GetObservedData();
        GetChiSquare(Chi_Square);
    }
    
    GetRotationCurve();

	return 0;
}
