#include "adapt_geometry.h"
#include "adapt_datastruct.h"
#include "adapt_compute.h"
#include <math.h>

double CheckFaceOrientation(std::vector<double> VcF, std::vector<std::vector<double> > Vfaces, std::vector<double> Vijk)
{
    std::vector<double> vnul(3);
    std::vector<double> vone(3);
    std::vector<double> rnul(3);
    
    double Lr0 = sqrt((VcF[0]-Vijk[0])*(VcF[0]-Vijk[0])
                      +(VcF[1]-Vijk[1])*(VcF[1]-Vijk[1])
                      +(VcF[2]-Vijk[2])*(VcF[2]-Vijk[2]));
    
    rnul[0] = (VcF[0]-Vijk[0])/Lr0;
    rnul[1] = (VcF[1]-Vijk[1])/Lr0;
    rnul[2] = (VcF[2]-Vijk[2])/Lr0;


    double orient0  = 1.0;
    int nppf = Vfaces.size();
    
    if(nppf==3) // triangle
    {
        vnul[0] = Vfaces[1][0]-Vfaces[0][0];
        vnul[1] = Vfaces[1][1]-Vfaces[0][1];
        vnul[2] = Vfaces[1][2]-Vfaces[0][2];

        vone[0] = Vfaces[2][0]-Vfaces[0][0];
        vone[1] = Vfaces[2][1]-Vfaces[0][1];
        vone[2] = Vfaces[2][2]-Vfaces[0][2];
        
        std::vector<double> n0 = ComputeSurfaceNormal(vnul,vone);
        orient0         = DotVec3D(rnul,n0);
    }
    if(nppf==4) // quad
    {
        vnul[0] = Vfaces[1][0]-Vfaces[0][0];
        vnul[1] = Vfaces[1][1]-Vfaces[0][1];
        vnul[2] = Vfaces[1][2]-Vfaces[0][2];

        vone[0] = Vfaces[3][0]-Vfaces[0][0];
        vone[1] = Vfaces[3][1]-Vfaces[0][1];
        vone[2] = Vfaces[3][2]-Vfaces[0][2];
        
        std::vector<double> n0 = ComputeSurfaceNormal(vnul,vone);
        orient0   = DotVec3D(rnul,n0);
    }
    
    return orient0;
}

Vec3D* ComputeSurfaceNormal_old(Vec3D* a, Vec3D* b)
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


std::vector<double> ComputeSurfaceNormal(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> V(3);
    std::vector<double> cross(3);

//        //Cross cross formula
//    cross[0] = (a->c1 * b->c2) - (a->c2 * b->c1);
//    cross[1] = (a->c2 * b->c0) - (a->c0 * b->c2);
//    cross[2] = (a->c0 * b->c1) - (a->c1 * b->c0);
    
    cross[0] = (a[1] * b[2]) - (a[2] * b[1]);
    cross[1] = (a[2] * b[0]) - (a[0] * b[2]);
    cross[2] = (a[0] * b[1]) - (a[1] * b[0]);

    double L = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    V[0]=cross[0]/L;
    V[1]=cross[1]/L;
    V[2]=cross[2]/L;
    
    return V;
}


double ComputeQuadSurfaceArea(double* P)
{
    
    std::vector<double> cross(3);
    
    std::vector<double> a(3);
    a[0] = P[1*3+0]-P[0*3+0];
    a[1] = P[1*3+1]-P[0*3+1];
    a[2] = P[1*3+2]-P[0*3+2];
    
    std::vector<double> b(3);
    b[0] = P[3*3+0]-P[0*3+0];
    b[1] = P[3*3+1]-P[0*3+1];
    b[2] = P[3*3+2]-P[0*3+2];
    
    //Cross cross formula
//    cross[0] = (a->c1 * b->c2) - (a->c2 * b->c1);
//    cross[1] = (a->c2 * b->c0) - (a->c0 * b->c2);
//    cross[2] = (a->c0 * b->c1) - (a->c1 * b->c0);
    
    cross[0] = (a[1] * b[2]) - (a[2] * b[1]);
    cross[1] = (a[2] * b[0]) - (a[0] * b[2]);
    cross[2] = (a[0] * b[1]) - (a[1] * b[0]);
    
    double R = fabs(sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]));

    return R;
}


double ComputeTriSurfaceArea(double* P)
{
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    std::vector<double> v2(3);
    
    v0[0] = P[0*3+0]; v1[0] = P[1*3+0];v2[0] = P[2*3+0];
    v0[1] = P[0*3+1]; v1[1] = P[1*3+1];v2[1] = P[2*3+1];
    v0[2] = P[0*3+2]; v1[2] = P[1*3+2];v2[2] = P[2*3+2];
    
    double a = ComputeEdgeLength(v0,v1);
    double b = ComputeEdgeLength(v0,v2);
    double c = ComputeEdgeLength(v1,v2);
    
    double p = (a+b+c)*0.5;
    
    double A = sqrt(p*(p-a)*(p-b)*(p-c));
    
    return A;
}

