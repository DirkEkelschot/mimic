#include "adapt_compute.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void NegateVec3D(Vec3D* a)
{
    a->c0 = -a->c0;
    a->c1 = -a->c1;
    a->c2 = -a->c2;
    
}

double DotVec3D(Vec3D* a, Vec3D* b)
{
    double res = a->c0*b->c0+a->c1*b->c1+a->c2*b->c2;
    
    return res;
}

Array<double>* MatMul(Array<double>* A, Array<double>* B)
{
    
    int n = A->getNrow();
    int o = A->getNcol();
    int m = B->getNcol();
    int k = B->getNrow();
    if(k!=o)
    {
        throw std::runtime_error("error :: Dimensions of A and B do not correspond.");
    }
    Array<double>* R = new Array<double>(n,m);
    
    double res = 0.0;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            res = 0.0;
            for(int k=0;k<o;k++)
            {
                res = res + A->getVal(i,k)*B->getVal(k,j);
            }
            R->setVal(i,j,res);
        }
    }
    return R;
}

double ComputeDetJac(double *P0,double *P1,double *P2,double *P3)
{
    
    double *JP1 = new double[9];
    
    JP1[0] = (P1[0]-P0[0]); JP1[1] = (P2[0]-P0[0]); JP1[2] = (P3[0]-P0[0]);
    JP1[3] = (P1[1]-P0[1]); JP1[4] = (P2[1]-P0[1]); JP1[5] = (P3[1]-P0[1]);
    JP1[6] = (P1[2]-P0[2]); JP1[7] = (P2[2]-P0[2]); JP1[8] = (P3[2]-P0[2]);
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
                 -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
                 +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    delete[] JP1;
    return DetJ;
}




double ComputeJ(double*P, int ElType){
    
    double J = 0.0;
    if (ElType==4)
    {
       double *P0  = new double[3];
       double *P1  = new double[3];
       double *P2  = new double[3];
       double *P3  = new double[3];
       
       P0[0]  = P[0*3+0]; P0[1]  = P[0*3+1]; P0[2] = P[0*3+2];
       P1[0]  = P[1*3+0]; P1[1]  = P[1*3+1]; P1[2] = P[1*3+2];
       P2[0]  = P[3*3+0]; P2[1]  = P[3*3+1]; P2[2] = P[3*3+2];
       P3[0]  = P[4*3+0]; P3[1]  = P[4*3+1]; P3[2] = P[4*3+2];
       double DJP0 = ComputeDetJac(P0,P1,P2,P3);

       P0[0]=P[1*3+0];P0[1]=P[1*3+1];P0[2]=P[1*3+2];
       P1[0]=P[2*3+0];P1[1]=P[2*3+1];P1[2]=P[2*3+2];
       P2[0]=P[0*3+0];P2[1]=P[0*3+1];P2[2]=P[0*3+2];
       P3[0]=P[5*3+0];P3[1]=P[5*3+1];P3[2]=P[5*3+2];
       double DJP1 = ComputeDetJac(P0,P1,P2,P3);
       
       // P2
       P0[0]=P[2*3+0];P0[1]=P[2*3+1];P0[2]=P[2*3+2];
       P1[0]=P[3*3+0];P1[1]=P[3*3+1];P1[2]=P[3*3+2];
       P2[0]=P[1*3+0];P2[1]=P[1*3+1];P2[2]=P[1*3+2];
       P3[0]=P[6*3+0];P3[1]=P[6*3+1];P3[2]=P[6*3+2];
       double DJP2 = ComputeDetJac(P0,P1,P2,P3);
       
       // P3
       P0[0]=P[3*3+0];P0[1]=P[3*3+1];P0[2]=P[3*3+2];
       P1[0]=P[0*3+0];P1[1]=P[0*3+1];P1[2]=P[0*3+2];
       P2[0]=P[2*3+0];P2[1]=P[2*3+1];P2[2]=P[2*3+2];
       P3[0]=P[7*3+0];P3[1]=P[7*3+1];P3[2]=P[7*3+2];
       double DJP3 = ComputeDetJac(P0,P1,P2,P3);
    
       // P4
       P0[0]=P[4*3+0];P0[1]=P[4*3+1];P0[2]=P[4*3+2];
       P1[0]=P[7*3+0];P1[1]=P[7*3+1];P1[2]=P[7*3+2];
       P2[0]=P[5*3+0];P2[1]=P[5*3+1];P2[2]=P[5*3+2];
       P3[0]=P[0*3+0];P3[1]=P[0*3+1];P3[2]=P[0*3+2];
       double DJP4 = ComputeDetJac(P0,P1,P2,P3);
       
       // P5
       P0[0]=P[5*3+0];P0[1]=P[5*3+1];P0[2]=P[5*3+2];
       P1[0]=P[4*3+0];P1[1]=P[4*3+1];P1[2]=P[4*3+2];
       P2[0]=P[6*3+0];P2[1]=P[6*3+1];P2[2]=P[6*3+2];
       P3[0]=P[1*3+0];P3[1]=P[1*3+1];P3[2]=P[1*3+2];
       double DJP5 = ComputeDetJac(P0,P1,P2,P3);
       
       // P6
       P0[0]=P[6*3+0];P0[1]=P[6*3+1];P0[2]=P[6*3+2];
       P1[0]=P[5*3+0];P1[1]=P[5*3+1];P1[2]=P[5*3+2];
       P2[0]=P[7*3+0];P2[1]=P[7*3+1];P2[2]=P[7*3+2];
       P3[0]=P[2*3+0];P3[1]=P[2*3+1];P3[2]=P[2*3+2];
       double DJP6 = ComputeDetJac(P0,P1,P2,P3);

       // P7
       P0[0]=P[7*3+0];P0[1]=P[7*3+1];P0[2]=P[7*3+2];
       P1[0]=P[6*3+0];P1[1]=P[6*3+1];P1[2]=P[6*3+2];
       P2[0]=P[4*3+0];P2[1]=P[4*3+1];P2[2]=P[4*3+2];
       P3[0]=P[3*3+0];P3[1]=P[3*3+1];P3[2]=P[3*3+2];
       double DJP7 = ComputeDetJac(P0,P1,P2,P3);
        
       delete[] P0;
       delete[] P1;
       delete[] P2;
       delete[] P3;
       J = (DJP0+DJP1+DJP2+DJP3+DJP4+DJP5+DJP6+DJP7)/8;

    }
    else if(ElType==6)
    {
        double *P0  = new double[3];
        double *P1  = new double[3];
        double *P2  = new double[3];
        double *P3  = new double[3];
        
        P0[0]  = P[0*3+0]; P0[1]  = P[0*3+1]; P0[2] = P[0*3+2];
        P1[0]  = P[1*3+0]; P1[1]  = P[1*3+1]; P1[2] = P[1*3+2];
        P2[0]  = P[2*3+0]; P2[1]  = P[2*3+1]; P2[2] = P[2*3+2];
        P3[0]  = P[3*3+0]; P3[1]  = P[3*3+1]; P3[2] = P[3*3+2];
        double DJP0 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[1*3+0]; P0[1]  = P[1*3+1]; P0[2] = P[1*3+2];
        P1[0]  = P[2*3+0]; P1[1]  = P[2*3+1]; P1[2] = P[2*3+2];
        P2[0]  = P[0*3+0]; P2[1]  = P[0*3+1]; P2[2] = P[0*3+2];
        P3[0]  = P[4*3+0]; P3[1]  = P[4*3+1]; P3[2] = P[4*3+2];
        double DJP1 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[2*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[0*3+0]; P1[1]  = P[0*3+1]; P1[2] = P[0*3+2];
        P2[0]  = P[1*3+0]; P2[1]  = P[1*3+1]; P2[2] = P[1*3+2];
        P3[0]  = P[5*3+0]; P3[1]  = P[5*3+1]; P3[2] = P[5*3+2];
        double DJP2 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[3*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[5*3+0]; P1[1]  = P[5*3+1]; P1[2] = P[5*3+2];
        P2[0]  = P[4*3+0]; P2[1]  = P[4*3+1]; P2[2] = P[4*3+2];
        P3[0]  = P[0*3+0]; P3[1]  = P[0*3+1]; P3[2] = P[0*3+2];
        double DJP3 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[4*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[3*3+0]; P1[1]  = P[3*3+1]; P1[2] = P[3*3+2];
        P2[0]  = P[5*3+0]; P2[1]  = P[5*3+1]; P2[2] = P[5*3+2];
        P3[0]  = P[1*3+0]; P3[1]  = P[1*3+1]; P3[2] = P[1*3+2];
        double DJP4 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[5*3+0]; P0[1]  = P[5*3+1]; P0[2] = P[5*3+2];
        P1[0]  = P[4*3+0]; P1[1]  = P[4*3+1]; P1[2] = P[4*3+2];
        P2[0]  = P[3*3+0]; P2[1]  = P[3*3+1]; P2[2] = P[3*3+2];
        P3[0]  = P[2*3+0]; P3[1]  = P[2*3+1]; P3[2] = P[2*3+2];
        double DJP5 = ComputeDetJac(P0,P1,P2,P3);
        
        delete[] P0;
        delete[] P1;
        delete[] P2;
        delete[] P3;
        
        J = (DJP0+DJP1+DJP2+DJP3+DJP4+DJP5)/6;
    }
    return J;
}


double ComputeEdgeLength(Vert* v0, Vert* v1)
{
    return sqrt((v0->x - v1->x) * (v0->x - v1->x)+
                (v0->y - v1->y) * (v0->y - v1->y)+
                (v0->z - v1->z) * (v0->z - v1->z));
}


double ComputeTetVolume(double *P)
{
    Vert* a = new Vert;
    Vert* b = new Vert;
    Vert* c = new Vert;
    Vert* d = new Vert;
    
    a->x = P[0*3+0]; b->x = P[1*3+0];
    a->y = P[0*3+1]; b->y = P[1*3+1];
    a->z = P[0*3+2]; b->z = P[1*3+2];
    
    c->x = P[2*3+0]; d->x = P[3*3+0];
    c->y = P[2*3+1]; d->y = P[3*3+1];
    c->z = P[2*3+2]; d->z = P[3*3+2];
    
    Vert* t = new Vert;
    t->x = (a->x-d->x);
    t->y = (a->y-d->y);
    t->z = (a->z-d->z);
    
    Vec3D* s = new Vec3D;
    s->c0 = (b->x-d->x);
    s->c1 = (b->y-d->y);
    s->c2 = (b->z-d->z);
    
    Vec3D* r= new Vec3D;
    r->c0 = (c->x-d->x);
    r->c1 = (c->y-d->y);
    r->c2 = (c->z-d->z);
    
    
    double cross[3];

        //Cross cross formula
    cross[0] = (s->c1 * r->c2) - (s->c2 * r->c1);
    cross[1] = (s->c2 * r->c0) - (s->c0 * r->c2);
    cross[2] = (s->c0 * r->c1) - (s->c1 * r->c0);
    
    double V = fabs(t->x*cross[0]+t->y*cross[1]+t->z*cross[2])/6.0;
    //delete t,s,r;
    delete t;
    delete s;
    delete r;
    delete a;
    delete b;
    delete c;
    delete d;
    //delete[] cross;
    return V;
    
    
}

double ComputeVolumeHexCell(double *P)
{
    
    double L01=0.0;
    double L15=0.0;
    double L04=0.0;
    double L45=0.0;
    double L37=0.0;
    double L23=0.0;
    double L26=0.0;
    double L67=0.0;
    
    double b0,b1,b2,b3;
    double H12=0.0,H47=0.0,H30=0.0,H56=0.0;
    
    Vert* v0 = new Vert;
    Vert* v1 = new Vert;
    
    v0->x = P[0*3+0]; v1->x = P[1*3+0];
    v0->y = P[0*3+1]; v1->y = P[1*3+1];
    v0->z = P[0*3+2]; v1->z = P[1*3+2];
    
    L01 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[1*3+0]; v1->x = P[5*3+0];
    v0->y = P[1*3+1]; v1->y = P[5*3+1];
    v0->z = P[1*3+2]; v1->z = P[5*3+2];
    
    L15 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[1*3+0]; v1->x = P[2*3+0];
    v0->y = P[1*3+1]; v1->y = P[2*3+1];
    v0->z = P[1*3+2]; v1->z = P[2*3+2];
    
    H12 = ComputeEdgeLength(v0,v1);
    b0 = 0.5*L01*L15;
    double vol0 = 1.0/3.0*b0*H12;
    //==================================================

    v0->x = P[0*3+0]; v1->x = P[4*3+0];
    v0->y = P[0*3+1]; v1->y = P[4*3+1];
    v0->z = P[0*3+2]; v1->z = P[4*3+2];
    
    L04 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[4*3+0]; v1->x = P[5*3+0];
    v0->y = P[4*3+1]; v1->y = P[5*3+1];
    v0->z = P[4*3+2]; v1->z = P[5*3+2];
    
    L45 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[4*3+0]; v1->x = P[7*3+0];
    v0->y = P[4*3+1]; v1->y = P[7*3+1];
    v0->z = P[4*3+2]; v1->z = P[7*3+2];
    
    H47 = ComputeEdgeLength(v0,v1);
    b1 = 0.5*L04*L45;
    double vol1 = 1.0/3.0*b1*H47;
    
    //==================================================
    
    v0->x = P[3*3+0]; v1->x = P[7*3+0];
    v0->y = P[3*3+1]; v1->y = P[7*3+1];
    v0->z = P[3*3+2]; v1->z = P[7*3+2];
    
    L37 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[2*3+0]; v1->x = P[3*3+0];
    v0->y = P[2*3+1]; v1->y = P[3*3+1];
    v0->z = P[2*3+2]; v1->z = P[3*3+2];
    
    L23 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[3*3+0]; v1->x = P[0*3+0];
    v0->y = P[3*3+1]; v1->y = P[0*3+1];
    v0->z = P[3*3+2]; v1->z = P[0*3+2];
   
    H30 = ComputeEdgeLength(v0,v1);
    b2 = 0.5*L37*L23;
    double vol2 = 1.0/3.0*b2*H30;
    
    //==================================================
    
    v0->x = P[2*3+0]; v1->x = P[6*3+0];
    v0->y = P[2*3+1]; v1->y = P[6*3+1];
    v0->z = P[2*3+2]; v1->z = P[6*3+2];
    
    L26 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[6*3+0]; v1->x = P[7*3+0];
    v0->y = P[6*3+1]; v1->y = P[7*3+1];
    v0->z = P[6*3+2]; v1->z = P[7*3+2];
    
    L67 = ComputeEdgeLength(v0,v1);
    
    v0->x = P[5*3+0]; v1->x = P[6*3+0];
    v0->y = P[5*3+1]; v1->y = P[6*3+1];
    v0->z = P[5*3+2]; v1->z = P[6*3+2];
    
    H56 = ComputeEdgeLength(v0,v1);
    b3 = 0.5*L26*L67;
    double vol3 = 1.0/3.0*b3*H56;
    delete v0;
    delete v1;
    return vol0+vol1+vol2+vol3;
}




double ComputeVolumeTetCell(double *P)
{
    Vert* a = new Vert;
    Vert* b = new Vert;
    Vert* c = new Vert;
    Vert* d = new Vert;
    
    a->x = P[0*3+0]; b->x = P[1*3+0];
    a->y = P[0*3+1]; b->y = P[1*3+1];
    a->z = P[0*3+2]; b->z = P[1*3+2];
    
    c->x = P[2*3+0]; d->x = P[3*3+0];
    c->y = P[2*3+1]; d->y = P[3*3+1];
    c->z = P[2*3+2]; d->z = P[3*3+2];
    
    Vert* t = new Vert;
    t->x = (a->x-d->x);
    t->y = (a->y-d->y);
    t->z = (a->z-d->z);
    
    Vec3D* s = new Vec3D;
    s->c0 = (b->x-d->x);
    s->c1 = (b->y-d->y);
    s->c2 = (b->z-d->z);
    
    Vec3D* r= new Vec3D;
    r->c0 = (c->x-d->x);
    r->c1 = (c->y-d->y);
    r->c2 = (c->z-d->z);
    
    
    double cross[3];

        //Cross cross formula
    cross[0] = (s->c1 * r->c2) - (s->c2 * r->c1);
    cross[1] = (s->c2 * r->c0) - (s->c0 * r->c2);
    cross[2] = (s->c0 * r->c1) - (s->c1 * r->c0);
    
    double V = fabs(t->x*cross[0]+t->y*cross[1]+t->z*cross[2])/6.0;
    
    delete a;
    delete b;
    delete c;
    delete d;
    delete t;
    delete s;
    delete r;
    
    return V;
}



double ComputeVolumePrismCell(double *P)
{
    Vert* v0 = new Vert;
    Vert* v1 = new Vert;
    Vert* v2 = new Vert;
    Vert* v3 = new Vert;
    Vert* v4 = new Vert;
    Vert* v5 = new Vert;
    
    v0->x = P[0*3+0]; v1->x = P[1*3+0];
    v0->y = P[0*3+1]; v1->y = P[1*3+1];
    v0->z = P[0*3+2]; v1->z = P[1*3+2];
    
    v2->x = P[2*3+0]; v3->x = P[3*3+0];
    v2->y = P[2*3+1]; v3->y = P[3*3+1];
    v2->z = P[2*3+2]; v3->z = P[3*3+2];
    
    v4->x = P[4*3+0]; v5->x = P[5*3+0];
    v4->y = P[4*3+1]; v5->y = P[5*3+1];
    v4->z = P[4*3+2]; v5->z = P[5*3+2];
    
    double La  = ComputeEdgeLength(v0,v2);
    
    double Lb  = ComputeEdgeLength(v1,v2);
    
    double Lc  = ComputeEdgeLength(v0,v1);
    
    double p  = La+Lb+Lc;
    //double pp = Lap+Lbp+Lcp;
    
    double A0  = sqrt(p/2.0*(p/2.0-La)*(p/2.0-Lb)*(p/2.0-Lc));
    //double A1 = sqrt(pp/2.0*(pp/2.0-Lap)*(pp/2.0-Lbp)*(pp/2.0-Lcp));
    
    double ha = ComputeEdgeLength(v1,v4);
    double hb = ComputeEdgeLength(v0,v3);
    double hc = ComputeEdgeLength(v2,v5);
    
//    double A2 = La*(hb+hc)/2.0;
//    double A3 = Lb*(ha+hc)/2.0;
//    double A4 = Lc*(ha+hb)/2.0;
//
//    double A = A0+A1+A2+A3+A4;
    double V = A0*(ha+hb+hc)/3.0;
    
    delete v0;
    delete v1;
    delete v2;
    delete v3;
    delete v4;
    delete v5;
    
    return V;
    
}





// This function outputs J as an array of 9 values where the matrix is defined as:

/*
 Jac = [J[0], J[1], J[2]
        J[3], J[4], J[5]
        J[6], J[7], J[8]]
*/
// J is computed using the 8-point isoparametric mapping for a hex. The 8-point rule should be sufficient since everything is linear anyways.
Vert* ComputeCenterCoord(double*P, int np)
{
    Vert* V = new Vert;
    
    int * ref = new int[np*3];
        
    // Allocate the arrays for the mapping function and its derivatives.
    double * N      = new double[np];
    
    for(int i=0;i<np;i++)
    {
        N[i]      = 0.0;
    }
    
    // Define the reference element as [-1,1]^3
    ref[0*3+0] = -1;    ref[0*3+1] = -1;    ref[0*3+2] = -1;
    ref[1*3+0] =  1;    ref[1*3+1] = -1;    ref[1*3+2] = -1;
    ref[2*3+0] =  1;    ref[2*3+1] =  1;    ref[2*3+2] = -1;
    ref[3*3+0] = -1;    ref[3*3+1] =  1;    ref[3*3+2] = -1;
            
    ref[4*3+0] = -1;    ref[4*3+1] = -1;    ref[4*3+2] = 1;
    ref[5*3+0] =  1;    ref[5*3+1] = -1;    ref[5*3+2] = 1;
    ref[6*3+0] =  1;    ref[6*3+1] =  1;    ref[6*3+2] = 1;
    ref[7*3+0] = -1;    ref[7*3+1] =  1;    ref[7*3+2] = 1;
     
    // We want to compute the Jacobian at the center of the cell
    // So we set eta, mu and ksi equal to 0 which is the center
    // of the reference cell.
    double eta = 0;
    double mu  = 0;
    double ksi = 0;
    V->x = 0.0;
    V->y = 0.0;
    V->z = 0.0;

    for(int i = 0; i < np; i++)
    {
        N[i] = 1.0/8.0*(1+ref[i*3+0]*eta)*(1+ref[i*3+1]*mu)*(1+ref[i*3+2]*ksi);
                        
        V->x = V->x+N[i]*P[i*3+0];
        V->y = V->y+N[i]*P[i*3+1];
        V->z = V->z+N[i]*P[i*3+2];
    }
    
    delete[] ref;
    delete[] N;
    return V;
}


Vert* ComputeCentroidCoord(double*P, int np)
{
    Vert* V = new Vert;
    //V = ComputeCenterCoord(P, np);
    V->x = 0.0;
    V->y = 0.0;
    V->z = 0.0;

    for(int i = 0; i < np; i++)
    {

        V->x = V->x+P[i*3+0];
        V->y = V->y+P[i*3+1];
        V->z = V->z+P[i*3+2];
    }

    V->x = V->x/np;
    V->y = V->y/np;
    V->z = V->z/np;
    
    return V;
}


double* ComputeJAtCenter(double*P, int np)
{
    //Vert V;
    double * J = new double[9];
    for(int i=0;i<9;i++)
    {
        J[i] = 0.0;
    }
    int * ref = new int[np*3];
        
    // Allocate the arrays for the mapping function and its derivatives.
    double * N      = new double[np];
    double * dNdeta = new double[np];
    double * dNdmu  = new double[np];
    double * dNdksi = new double[np];
    
    for(int i=0;i<np;i++)
    {
        N[i]      = 0.0;
        dNdeta[i] = 0.0;
        dNdmu[i]  = 0.0;
        dNdksi[i] = 0.0;
    }
    
    // Define the reference element as [-1,1]^3
    ref[0*3+0] = -1;    ref[0*3+1] = -1;    ref[0*3+2] = -1;
    ref[1*3+0] =  1;    ref[1*3+1] = -1;    ref[1*3+2] = -1;
    ref[2*3+0] =  1;    ref[2*3+1] =  1;    ref[2*3+2] = -1;
    ref[3*3+0] = -1;    ref[3*3+1] =  1;    ref[3*3+2] = -1;
            
    ref[4*3+0] = -1;    ref[4*3+1] = -1;    ref[4*3+2] = 1;
    ref[5*3+0] =  1;    ref[5*3+1] = -1;    ref[5*3+2] = 1;
    ref[6*3+0] =  1;    ref[6*3+1] =  1;    ref[6*3+2] = 1;
    ref[7*3+0] = -1;    ref[7*3+1] =  1;    ref[7*3+2] = 1;
     
    // We want to compute the Jacobian at the center of the cell
    // So we set eta, mu and ksi equal to 0 which is the center
    // of the reference cell.
    double eta = 0;
    double mu  = 0;
    double ksi = 0;
    double xphys = 0.0;double yphys = 0.0;double zphys = 0.0;

    // The basis functions for the 8 point isoparametric mapping for a hex
    // looks like:
    
    /*
    N[0] = 1.0/8.0*(1-eta)*(1-mu)*(1-ksi);
    N[1] = 1.0/8.0*(1+eta)*(1-mu)*(1-ksi);
    N[2] = 1.0/8.0*(1+eta)*(1+mu)*(1-ksi);
    N[3] = 1.0/8.0*(1-eta)*(1+mu)*(1-ksi);
    N[4] = 1.0/8.0*(1-eta)*(1-mu)*(1+ksi);
    N[5] = 1.0/8.0*(1+eta)*(1-mu)*(1+ksi);
    N[6] = 1.0/8.0*(1+eta)*(1+mu)*(1+ksi);
    N[7] = 1.0/8.0*(1-eta)*(1+mu)*(1+ksi);
    */
        
    for(int i = 0; i < np; i++)
    {
        N[i] = 1.0/8.0*(1+ref[i*3+0]*eta)*(1+ref[i*3+1]*mu)*(1+ref[i*3+2]*ksi);
        
        dNdeta[i] = (1.0/8.0 * (1+ref[i*3+1]*mu)  * (1+ref[i*3+2]*ksi))*ref[i*3+0];
        dNdmu[i]  = (1.0/8.0 * (1+ref[i*3+0]*eta) * (1+ref[i*3+2]*ksi))*ref[i*3+1];
        dNdksi[i] = (1.0/8.0 * (1+ref[i*3+0]*eta) * (1+ref[i*3+1]*mu))*ref[i*3+2];
                        
        xphys = xphys+N[i]*P[i*3+0];
        yphys = yphys+N[i]*P[i*3+1];
        zphys = zphys+N[i]*P[i*3+2];
            
        J[0] = J[0]+dNdeta[i]*P[i*3+0];
        J[1] = J[1]+dNdeta[i]*P[i*3+1];
        J[2] = J[2]+dNdeta[i]*P[i*3+2];
        
        J[3] = J[3]+dNdmu[i]*P[i*3+0];
        J[4] = J[4]+dNdmu[i]*P[i*3+1];
        J[5] = J[5]+dNdmu[i]*P[i*3+2];
            
        J[6] = J[6]+dNdksi[i]*P[i*3+0];
        J[7] = J[7]+dNdksi[i]*P[i*3+1];
        J[8] = J[8]+dNdksi[i]*P[i*3+2];
            
    }
    
    return J;
}

//double* ComputeJAtCenter_tet(double*P)
//{
//    Array<double>* ref = new Array<double>(4,4);
//    ref->setVal(0,0,-1.0);ref->setVal(0,1, 1.0);ref->setVal(0,2,-1.0);ref->setVal(0,3,-1.0);
//    ref->setVal(1,0,-1.0);ref->setVal(1,1,-1.0);ref->setVal(1,2, 1.0);ref->setVal(1,3,-1.0);
//    ref->setVal(2,0,-1.0);ref->setVal(2,1,-1.0);ref->setVal(2,2,-1.0);ref->setVal(2,3, 1.0);
//    ref->setVal(3,0, 1.0);ref->setVal(3,1, 1.0);ref->setVal(3,2, 1.0);ref->setVal(3,3, 1.0);
//
//    Array<double>* phys = new Array<double>(4,4);
//    phys->setVal(0,0,P[0*3+0]);phys->setVal(0,1,P[1*3+0]);phys->setVal(0,2,P[2*3+0]);phys->setVal(0,3,P[3*3+0]);
//    phys->setVal(1,0,P[0*3+1]);phys->setVal(1,1,P[1*3+1]);phys->setVal(1,2,P[2*3+1]);phys->setVal(1,3,P[3*3+1]);
//    phys->setVal(2,0,P[0*3+2]);phys->setVal(2,1,P[1*3+2]);phys->setVal(2,2,P[2*3+2]);phys->setVal(2,3,P[3*3+2]);
//    phys->setVal(3,0,    1.0); phys->setVal(3,1,     1.0);phys->setVal(3,2,     1.0);phys->setVal(3,3,     1.0);
//
//    Array<double>* phys_t = new Array<double>(4,4);
//    Array<double>* ref_t = new Array<double>(4,4);
//
//    for(int i=0;i<4;i++)
//    {
//        for(int j=0;j<4;j++)
//        {
//            ref_t->setVal(j,i,ref->getVal(i,j));
//            phys_t->setVal(j,i,phys->getVal(i,j));
//        }
//    }
//
//    Array<double>* iphys = MatInv(phys);
//    Array<double>* iref  = MatInv(ref);
//    Array<double>* M  = MatMul(ref,iphys);
//    Array<double>* Ma = MatMul(phys,iref);
//
//    double* J = new double[3*3];
//    for(int i=0;i<3;i++)
//    {
//        for(int j=0;j<3;j++)
//        {
//            J[i*3+j] = Ma->getVal(i,j);
//        }
//    }
//
//    return J;
//}


Array<double>* ComputeJAtCenter_tet_v2(double*P)
{

    Array<double>* a = new Array<double>(3,3);
    a->setVal(0,0,P[1*3+0]-P[0*3+0]);
    a->setVal(0,1,P[2*3+0]-P[0*3+0]);
    a->setVal(0,2,P[3*3+0]-P[0*3+0]);
    
    a->setVal(1,0,P[1*3+0]-P[0*3+1]);
    a->setVal(1,1,P[2*3+0]-P[0*3+1]);
    a->setVal(1,2,P[3*3+0]-P[0*3+1]);
    
    a->setVal(2,0,P[1*3+0]-P[0*3+2]);
    a->setVal(2,1,P[2*3+0]-P[0*3+2]);
    a->setVal(2,2,P[3*3+0]-P[0*3+2]);
    
    Array<double>* w0 = new Array<double>(3,3);
    w0->setVal(0,0,1.0);
    w0->setVal(0,1,0.5);
    w0->setVal(0,2,0.5);
    
    w0->setVal(1,0,0.0);
    w0->setVal(1,1,sqrt(3.0)/2.0);
    w0->setVal(1,2,sqrt(3.0)/6.0);
    
    w0->setVal(2,0,0.0);
    w0->setVal(2,1,0.0);
    w0->setVal(2,2,sqrt(2.0)/sqrt(3.0));

    Array<double>* iw0 = MatInv(w0);
    Array<double>* M  = MatMul(a,iw0);
    delete a;
    delete iw0;
    
    return M;
}
double ComputeDeterminantJ_tet_v2(double*P)
{
    Array<double>* JP1 = ComputeJAtCenter_tet_v2(P);

    double J0 = JP1->getVal(0,0);
    double J1 = JP1->getVal(1,0);
    double J2 = JP1->getVal(2,0);
    double J3 = JP1->getVal(0,1);
    double J4 = JP1->getVal(1,1);
    double J5 = JP1->getVal(2,1);
    double J6 = JP1->getVal(0,2);
    double J7 = JP1->getVal(1,2);
    double J8 = JP1->getVal(2,2);
    
    double DetJ = J0*(J4*J8-J7*J5)
    -J1*(J3*J8-J6*J5)
    +J2*(J3*J7-J6*J4);
    
    delete JP1;
    
    return DetJ;
}

double ComputeDeterminantJ_tet(double*P)
{
    double* J0=new double[9];
    std::vector<double> dJvec;
    J0[0] = P[1*3+0]-P[0*3+0];
    J0[1] = P[2*3+0]-P[0*3+0];
    J0[2] = P[3*3+0]-P[0*3+0];

    J0[3] = P[1*3+1]-P[0*3+1];
    J0[4] = P[2*3+1]-P[0*3+1];
    J0[5] = P[3*3+1]-P[0*3+1];

    J0[6] = P[1*3+2]-P[0*3+2];
    J0[7] = P[2*3+2]-P[0*3+2];
    J0[8] = P[3*3+2]-P[0*3+2];
    
    double DetJ0 =   J0[0]*(J0[4]*J0[8]-J0[7]*J0[5])
                    -J0[1]*(J0[3]*J0[8]-J0[6]*J0[5])
                    +J0[2]*(J0[3]*J0[7]-J0[6]*J0[4]);
    
    dJvec.push_back(DetJ0);

    J0[0] = P[0*3+0]-P[1*3+0];
    J0[1] = P[2*3+0]-P[1*3+0];
    J0[2] = P[3*3+0]-P[1*3+0];

    J0[3] = P[0*3+1]-P[1*3+1];
    J0[4] = P[2*3+1]-P[1*3+1];
    J0[5] = P[3*3+1]-P[1*3+1];

    J0[6] = P[0*3+2]-P[1*3+2];
    J0[7] = P[2*3+2]-P[1*3+2];
    J0[8] = P[3*3+2]-P[1*3+2];
    
    DetJ0 =  J0[0]*(J0[4]*J0[8]-J0[7]*J0[5])
            -J0[1]*(J0[3]*J0[8]-J0[6]*J0[5])
            +J0[2]*(J0[3]*J0[7]-J0[6]*J0[4]);
    
    dJvec.push_back(-DetJ0);
    
    J0[0] = P[0*3+0]-P[2*3+0];
    J0[1] = P[1*3+0]-P[2*3+0];
    J0[2] = P[3*3+0]-P[2*3+0];

    J0[3] = P[0*3+1]-P[2*3+1];
    J0[4] = P[1*3+1]-P[2*3+1];
    J0[5] = P[3*3+1]-P[2*3+1];

    J0[6] = P[0*3+2]-P[2*3+2];
    J0[7] = P[1*3+2]-P[2*3+2];
    J0[8] = P[3*3+2]-P[2*3+2];
    
    DetJ0 =  J0[0]*(J0[4]*J0[8]-J0[7]*J0[5])
            -J0[1]*(J0[3]*J0[8]-J0[6]*J0[5])
            +J0[2]*(J0[3]*J0[7]-J0[6]*J0[4]);
    
    dJvec.push_back(DetJ0);
    
    J0[0] = P[0*3+0]-P[3*3+0];
    J0[1] = P[1*3+0]-P[3*3+0];
    J0[2] = P[2*3+0]-P[3*3+0];

    J0[3] = P[0*3+1]-P[3*3+1];
    J0[4] = P[1*3+1]-P[3*3+1];
    J0[5] = P[2*3+1]-P[3*3+1];

    J0[6] = P[0*3+2]-P[3*3+2];
    J0[7] = P[1*3+2]-P[3*3+2];
    J0[8] = P[2*3+2]-P[3*3+2];
    
    DetJ0 =  J0[0]*(J0[4]*J0[8]-J0[7]*J0[5])
            -J0[1]*(J0[3]*J0[8]-J0[6]*J0[5])
            +J0[2]*(J0[3]*J0[7]-J0[6]*J0[4]);
    
    dJvec.push_back(-DetJ0);
    //std::cout << dJvec[0] << " " <<  dJvec[1] << " " << dJvec[2] << " " << dJvec[3] << std::endl;
    double DetJ = *std::min_element(dJvec.begin(),dJvec.end());
    delete[] J0;
    return DetJ;
}


double ComputeDeterminantJ(double*P, int np)
{
    double* JP1 = ComputeJAtCenter(P, np);
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
    -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
    +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    return DetJ;
    
}


Array<double>* ComputeHessian(Partition_old* pa)
{
    int ndim = pa->ndim;
    Array<double>* H = new Array<double>(pa->ien->getNrow(),ndim);
    
    return H;
}


Array<double>* ComputeDeterminantofJacobian(ParArray<int>* ien, Array<double>* xcn)
{
    
    int np = 8; // Assuming its all hexes.
    int nel_loc = ien->getNrow();
    int ncol = ien->getNcol();
    Array<double>* detJ = new Array<double>(nel_loc,1);
    double* P = new double[np*3];
    int Vid, Vlid;
    for(int i=0;i<nel_loc;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            Vid = ien->getVal(i,j);
            
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        double dJ = ComputeDeterminantJ(P,8);
        
        detJ->setVal(i,0,dJ);
        
        //std::cout << "================================" << std::endl;
        //std::cout << J[0] << " " << J[1] << " " << J[2] << std::endl;
        //std::cout << J[3] << " " << J[4] << " " << J[5] << std::endl;
        //std::cout << J[6] << " " << J[7] << " " << J[8] << std::endl;
        //std::cout << "================================" << std::endl;
    }
    delete[] P;
    
    return detJ;
}


double* ComputeVolumeCellsSerial(Array<double>* xcn, Array<int>* ien, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int Nel = ien->getNrow();
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(Nel/world_size) + MIN(world_rank, Nel%world_size);
    
    int Nelements = Nel;
    double * vol_cells = new double[Nelements];
    int np = 8;
    double* P = new double[np*3];
    
    
    double vhex = 0.0;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i+offset,j+1)-1;
            
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        vhex = ComputeVolumeHexCell(P);
    
        vol_cells[i] = vhex;
        
    }
    
    delete[] P;
    
    return vol_cells;
    
}







double* ComputeVolumeCellsReducedToVerts(Array<double>* xcn, Array<int>* ien)
{
    
    
    int Nelements = ien->getNrow();
    int Nnodes = xcn->getNrow();
    
    double * vert_cnt   = new double[Nnodes];
    double * sum_vol    = new double[Nnodes];
    double * vol_verts  = new double[Nnodes];
    
    
    int np = 8;
    double* P = new double[np*3];
    
    double vhex = 0.0;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid      = ien->getVal(i,j+1)-1;
            
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        vhex = ComputeVolumeHexCell(P);
        
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i,j+1)-1;
            
            vert_cnt[Vid] = vert_cnt[Vid] + 1;
            sum_vol[Vid]  = sum_vol[Vid] + vhex;
        }
    }
    
    for(int i=0;i<Nnodes;i++)
    {
        vol_verts[i] = sum_vol[i]/vert_cnt[i];
    }
    
    
    delete[] vert_cnt;
    delete[] sum_vol;
    delete[] P;
    
    return vol_verts;
}



void UnitTestJacobian()
{
    double* Hex = new double[8*3];
    
    Hex[0*3+0] = 0;     Hex[0*3+1] = 0;     Hex[0*3+2] = 0;
    Hex[1*3+0] = 0.5;   Hex[1*3+1] = 0;     Hex[1*3+2] = 0;
    Hex[2*3+0] = 0.5;   Hex[2*3+1] = 0.5;   Hex[2*3+2] = 0;
    Hex[3*3+0] = 0;     Hex[3*3+1] = 0.5;   Hex[3*3+2] = 0;
    
    Hex[4*3+0] = 0;     Hex[4*3+1] = 0;     Hex[4*3+2] = 0.5;
    Hex[5*3+0] = 0.5;   Hex[5*3+1] = 0;     Hex[5*3+2] = 0.5;
    Hex[6*3+0] = 0.5;   Hex[6*3+1] = 0.5;   Hex[6*3+2] = 0.5;
    Hex[7*3+0] = 0;     Hex[7*3+1] = 0.5;   Hex[7*3+2] = 0.5;
    
    double* JP1 = ComputeJAtCenter(Hex,8);
    
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << JP1[i*3+j] << " ";
        }
        std::cout << " " << std::endl;
    }
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
                 -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
                 +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    std::cout << DetJ << std::endl;
}







void ComputeMetric(Partition* Pa,
                        std::vector<double> metric_inputs, MPI_Comm comm, std::map<int,Array<double>* > &Hess_vm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++
    double hmin         = metric_inputs[1];
    double hmax         = metric_inputs[2];
    double f            = metric_inputs[3];
    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double* Hmet = new double[9];
    

    std::map<int,Array<double>*>::iterator itm;
    int i = 0;
    for(itm=Hess_vm.begin();itm!=Hess_vm.end();itm++)
    {
        int glob_vid = itm->first;
        
        Hmet[0] = Hess_vm[glob_vid]->getVal(0,0);
        Hmet[1] = Hess_vm[glob_vid]->getVal(1,0);
        Hmet[2] = Hess_vm[glob_vid]->getVal(2,0);

        Hmet[3] = Hess_vm[glob_vid]->getVal(1,0);
        Hmet[4] = Hess_vm[glob_vid]->getVal(3,0);
        Hmet[5] = Hess_vm[glob_vid]->getVal(4,0);

        Hmet[6] = Hess_vm[glob_vid]->getVal(2,0);
        Hmet[7] = Hess_vm[glob_vid]->getVal(4,0);
        Hmet[8] = Hess_vm[glob_vid]->getVal(5,0);
        
        //delete Hess_vm[glob_vid];
        
        double * WR = new double[3];
        double * WI = new double[3];
        double * V  = new double[3*3];
        double * iV = new double[3*3];
        double* WRn = new double[3];

        Array<double>* DR  = new Array<double>(3,3);
        Array<double>* UR  = new Array<double>(3,3);
        for(int j=0;j<3;j++)
        {
            WR[j]  = 0.0;
            WRn[j] = 0.0;
            WI[j]  = 0.0;

            for(int k=0;k<3;k++)
            {
                DR->setVal(j,k,0.0);
                V[j*3+k] = 0.0;
                iV[j*3+k] = 0.0;
            }
        }
        
        Eig* eig = ComputeEigenDecomp(3, Hmet);
        
        WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[2] = std::min(std::max(f*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        
        Hess_vm[glob_vid] = Rf;
        
        delete DR;
        delete UR;
        delete Rs;
        //delete Rf;

        delete[] iV;
        delete[] V;
        delete[] WR;
        delete[] WI;
        delete[] WRn;
        i++;
    }
    
    delete[] Hmet;
}


Array<double>* ComputeFaceValues(Partition* P, Array<double>* U, MPI_Comm comm)
{
    int nface, start, end, rank, size;
    int nloc, offset, adjEl_id, leid, geid, i, t;
    double u_o,u_c;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    //offset  = P->getLocalPartition()->getOffset(rank);
    int* xadj   = P->getXadj();
    int* adjcny = P->getAdjcny();
    nloc    = P->getLocalPartition()->getNrow();
    
    nface = 6; // # hardcoded for hexes for now
    
    std::vector<double> Uelem_all         = P->PartitionAuxilaryData(U, comm);
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::map<int,int> lE2gE               = P->getLocalElement2GlobalElement();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    Array<double>* Uf                     = new Array<double>(nloc,nface);
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();

    
    
    for(i=0;i<nloc;i++)
    {
        start  = xadj[i];
        end    = xadj[i+1];
        u_c    = U->getVal(i,0);
        t = 0;
        for(int j=start;j<end;j++)
        {
            adjEl_id = adjcny[j];
            leid     = gE2lE[adjEl_id];
            u_o      = Uelem_all[leid];
            
            Uf->setVal(i,t,u_c-u_o);
            
            t++;
        }
    }
    
    return Uf;
}


Array<double>* ComputeVolumes(Partition* Pa)
{
    int i,j,k;
    int loc_vid;
    std::map<int,int> gE2lE                 = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem               = Pa->getLocElem();
    int nLocElem                            = Loc_Elem.size();
    std::map<int,std::vector<int> > gE2lV   = Pa->getGlobElem2LocVerts();
    std::vector<Vert*> locVerts              = Pa->getLocalVerts();
    Array<double>* Volumes = new Array<double>(nLocElem,1);
    std::vector<int> vijkIDs;
    double* Pijk = new double[8*3];
    for(i=0;i<nLocElem;i++)
    {
       int gEl = Loc_Elem[i];

       vijkIDs = gE2lV[gEl];

       for(k=0;k<vijkIDs.size();k++)
       {
          loc_vid     = vijkIDs[k];
          Pijk[k*3+0] = locVerts[loc_vid]->x;
          Pijk[k*3+1] = locVerts[loc_vid]->y;
          Pijk[k*3+2] = locVerts[loc_vid]->z;
       }

       double Vol = ComputeVolumeHexCell(Pijk);
       Volumes->setVal(i,0,Vol);
    }
    
    delete[] Pijk;
    return Volumes;
}


