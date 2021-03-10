#include "adapt_geometry.h"
#include "adapt_datatype.h"
#include "adapt_datastruct.h"
#include "adapt_compute.h"
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


double ComputeQuadSurfaceArea(double* P)
{
    
    double cross[3];
    Vec3D* a = new Vec3D;
    a->c0 = P[1*3+0]-P[0*3+0];
    a->c1 = P[1*3+1]-P[0*3+1];
    a->c2 = P[1*3+2]-P[0*3+2];
    Vec3D* b = new Vec3D;
    b->c0 = P[3*3+0]-P[0*3+0];
    b->c1 = P[3*3+1]-P[0*3+1];
    b->c2 = P[3*3+2]-P[0*3+2];
    //Cross cross formula
    cross[0] = (a->c1 * b->c2) - (a->c2 * b->c1);
    cross[1] = (a->c2 * b->c0) - (a->c0 * b->c2);
    cross[2] = (a->c0 * b->c1) - (a->c1 * b->c0);
    
    double R = fabs(sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));
    delete a;
    delete b;
    return R;
}


double ComputeTriSurfaceArea(double* P)
{

    Vert* v0 = new Vert;
    Vert* v1 = new Vert;
    Vert* v2 = new Vert;
    
    v0->x = P[0*3+0]; v1->x = P[1*3+0];v2->x = P[2*3+0];
    v0->y = P[0*3+1]; v1->y = P[1*3+1];v2->y = P[2*3+1];
    v0->z = P[0*3+2]; v1->z = P[1*3+2];v2->z = P[2*3+2];
    
    double a = ComputeEdgeLength(v0,v1);
    double b = ComputeEdgeLength(v0,v2);
    double c = ComputeEdgeLength(v1,v2);
    
    double p = (a+b+c)*0.5;
    
    double A = sqrt(p*(p-a)*(p-b)*(p-c));
    
    delete v0,v1,v2;
    
    return A;
}

