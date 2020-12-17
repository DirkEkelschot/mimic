#include "adapt_geometry.h"
#include <math.h>


Vec3D* ComputeSurfaceNormal(Vec3D* a, Vec3D* b)
{
    Vec3D* V = new Vec3D;
    double cross[3];

        //Cross cross formula
    cross[0] = (a->c1 * b->c2) - (a->c2 * b->c1);
    cross[1] = (a->c2 * b->c0) - (a->c0 * b->c2);
    cross[2] = (a->c0 * b->c1) - (a->c1 * b->c0);

    double L = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    V->c0=cross[0]/L;
    V->c1=cross[1]/L;
    V->c2=cross[2]/L;
    return V;
}


double ComputeSurfaceArea(Vec3D* a, Vec3D* b)
{
    double cross[3];

    //Cross cross formula
    cross[0] = (a->c1 * b->c2) - (a->c2 * b->c1);
    cross[1] = (a->c2 * b->c0) - (a->c0 * b->c2);
    cross[2] = (a->c0 * b->c1) - (a->c1 * b->c0);
    
    double R = fabs(sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
    return R;
}

