#include "adapt_geometry.h"
#include <math.h>


Vec3D* ComputeSurfaceNormal(double *a, double *b)
{
    Vec3D* V = new Vec3D;
    double cross[3];

        //Cross cross formula
    cross[0] = (a[1] * b[2]) - (a[2] * b[1]);
    cross[1] = (a[2] * b[0]) - (a[0] * b[2]);
    cross[2] = (a[0] * b[1]) - (a[1] * b[0]);

    double L = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    V->c0=cross[0]/L;
    V->c1=cross[1]/L;
    V->c2=cross[2]/L;
    
    return V;
}


double ComputeSurfaceArea(double *a, double *b)
{
    double cross[3];

    //Cross cross formula
    cross[0] = (a[1] * b[2]) - (a[2] * b[1]);
    cross[1] = (a[2] * b[0]) - (a[0] * b[2]);
    cross[2] = (a[0] * b[1]) - (a[1] * b[0]);
    
    return fabs(sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
}

