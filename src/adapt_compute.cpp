#include "adapt_compute.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


std::vector<std::vector<double> > MatInv_Lite(std::vector<std::vector<double> > A)
{
    int n = A.size();
    int size = n*n;
    double WORK [size];
    int info;
    int Pivot[n];
    std::vector<double> R_tmp(n*n,0.0);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            R_tmp[n*i+j] = A[i][j];
        }
    }
    
    dgetrf_(&n, &n, R_tmp.data(), &n, Pivot, &info);
    dgetri_(&n, R_tmp.data(), &n, Pivot, WORK, &size, &info);

    std::vector<std::vector<double> > R(n);

    for(int i=0;i<n;i++)
    {
        std::vector<double> row(n,0.0);
        for(int j=0;j<3;j++)
        {
            row[j] = R_tmp[i*n+j];
        }

        R[i] = row;
    }
    
    return R;
}

std::vector<std::vector<double> > MatMul_Lite(std::vector<std::vector<double> > A, 
                           std::vector<std::vector<double> > B)
{
    
    int n = A.size();
    int o = A[0].size();
    int m = B.size();
    int k = B[0].size();

    if(k!=o)
    {
        throw std::runtime_error("error :: Dimensions of A and B do not correspond.");
    }

    std::vector<std::vector<double> > R(n);
    //Array<double>* R = new Array<double>(n,m);
    
    double res = 0.0;
    for(int i=0;i<n;i++)
    {
        std::vector<double> row(m,0.0);

        for(int j=0;j<m;j++)
        {
            res = 0.0;
            for(int k=0;k<o;k++)
            {
                res = res + A[i][k]*B[k][j];
            }
            row[j] = res;
        }

        R[i] = row;

    }
    return R;
}



void NegateVec3D(std::vector<double> a)
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
    
}

double DotVec3D(std::vector<double> a, std::vector<double> b)
{
    double res = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
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


double ComputeEdgeLength(std::vector<double> v0, std::vector<double> v1)
{
    return sqrt((v0[0] - v1[0]) * (v0[0] - v1[0])+
                (v0[1] - v1[1]) * (v0[1] - v1[1])+
                (v0[2] - v1[2]) * (v0[2] - v1[2]));
}


double ComputeTetVolume(std::vector<double> P)
{
    std::vector<double> a(3);
    std::vector<double> b(3);
    std::vector<double> c(3);
    std::vector<double> d(3);
    
    a[0] = P[0*3+0]; b[0] = P[1*3+0];
    a[1] = P[0*3+1]; b[1] = P[1*3+1];
    a[2] = P[0*3+2]; b[2] = P[1*3+2];
    
    c[0] = P[2*3+0]; d[0] = P[3*3+0];
    c[1] = P[2*3+1]; d[1] = P[3*3+1];
    c[2] = P[2*3+2]; d[2] = P[3*3+2];
    
    
    std::vector<double> t(3);
    t[0] = (a[0]-d[0]);
    t[1] = (a[1]-d[1]);
    t[2] = (a[2]-d[2]);
    
    
    std::vector<double> s(3);
    s[0] = (b[0]-d[0]);
    s[1] = (b[1]-d[1]);
    s[2] = (b[2]-d[2]);
    
    
    std::vector<double> r(3);
    r[0] = (c[0]-d[0]);
    r[1] = (c[1]-d[1]);
    r[2] = (c[2]-d[2]);
    
    
    double cross[3];

        //Cross cross formula
    cross[0] = (s[1] * r[2]) - (s[2] * r[1]);
    cross[1] = (s[2] * r[0]) - (s[0] * r[2]);
    cross[2] = (s[0] * r[1]) - (s[1] * r[0]);
    
    double V = fabs(t[0]*cross[0]+t[1]*cross[1]+t[2]*cross[2])/6.0;
    //delete t,s,r;

    //delete[] cross;
    return V;
    
    
}

double ComputeVolumeHexCell(std::vector<double> P)
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
    
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    
    v0[0] = P[0*3+0]; v1[0] = P[1*3+0];
    v0[1] = P[0*3+1]; v1[1] = P[1*3+1];
    v0[2] = P[0*3+2]; v1[2] = P[1*3+2];
    
    L01 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[1*3+0]; v1[0] = P[5*3+0];
    v0[1] = P[1*3+1]; v1[1] = P[5*3+1];
    v0[2] = P[1*3+2]; v1[2] = P[5*3+2];
    
    L15 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[1*3+0]; v1[0] = P[2*3+0];
    v0[1] = P[1*3+1]; v1[1] = P[2*3+1];
    v0[2] = P[1*3+2]; v1[2] = P[2*3+2];
    
    H12 = ComputeEdgeLength(v0,v1);
    b0 = 0.5*L01*L15;
    double vol0 = 1.0/3.0*b0*H12;
    //==================================================

    v0[0] = P[0*3+0]; v1[0] = P[4*3+0];
    v0[1] = P[0*3+1]; v1[1] = P[4*3+1];
    v0[2] = P[0*3+2]; v1[2] = P[4*3+2];
    
    L04 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[4*3+0]; v1[0] = P[5*3+0];
    v0[1] = P[4*3+1]; v1[1] = P[5*3+1];
    v0[2] = P[4*3+2]; v1[2] = P[5*3+2];
    
    L45 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[4*3+0]; v1[0] = P[7*3+0];
    v0[1] = P[4*3+1]; v1[1] = P[7*3+1];
    v0[2] = P[4*3+2]; v1[2] = P[7*3+2];
    
    H47 = ComputeEdgeLength(v0,v1);
    b1 = 0.5*L04*L45;
    double vol1 = 1.0/3.0*b1*H47;
    
    //==================================================
    
    v0[0] = P[3*3+0]; v1[0] = P[7*3+0];
    v0[1] = P[3*3+1]; v1[1] = P[7*3+1];
    v0[2] = P[3*3+2]; v1[2] = P[7*3+2];
    
    L37 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[2*3+0]; v1[0] = P[3*3+0];
    v0[1] = P[2*3+1]; v1[1] = P[3*3+1];
    v0[2] = P[2*3+2]; v1[2] = P[3*3+2];
    
    L23 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[3*3+0]; v1[0] = P[0*3+0];
    v0[1] = P[3*3+1]; v1[1] = P[0*3+1];
    v0[2] = P[3*3+2]; v1[2] = P[0*3+2];
   
    H30 = ComputeEdgeLength(v0,v1);
    b2 = 0.5*L37*L23;
    double vol2 = 1.0/3.0*b2*H30;
    
    //==================================================
    
    v0[0] = P[2*3+0]; v1[0] = P[6*3+0];
    v0[1] = P[2*3+1]; v1[1] = P[6*3+1];
    v0[2] = P[2*3+2]; v1[2] = P[6*3+2];
    
    L26 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[6*3+0]; v1[0] = P[7*3+0];
    v0[1] = P[6*3+1]; v1[1] = P[7*3+1];
    v0[2] = P[6*3+2]; v1[2] = P[7*3+2];
    
    L67 = ComputeEdgeLength(v0,v1);
    
    v0[0] = P[5*3+0]; v1[0] = P[6*3+0];
    v0[1] = P[5*3+1]; v1[1] = P[6*3+1];
    v0[2] = P[5*3+2]; v1[2] = P[6*3+2];
    
    H56 = ComputeEdgeLength(v0,v1);
    b3 = 0.5*L26*L67;
    double vol3 = 1.0/3.0*b3*H56;

    return vol0+vol1+vol2+vol3;
}




double ComputeVolumeTetCell(std::vector<double> P)
{

    
    std::vector<double> a(3);
    std::vector<double> b(3);
    std::vector<double> c(3);
    std::vector<double> d(3);
    
    a[0] = P[0*3+0]; b[0] = P[1*3+0];
    a[1] = P[0*3+1]; b[1] = P[1*3+1];
    a[2] = P[0*3+2]; b[2] = P[1*3+2];
    
    c[0] = P[2*3+0]; d[0] = P[3*3+0];
    c[1] = P[2*3+1]; d[1] = P[3*3+1];
    c[2] = P[2*3+2]; d[2] = P[3*3+2];
    
    std::vector<double> t(3);
    t[0] = (a[0]-d[0]);
    t[1] = (a[1]-d[1]);
    t[2] = (a[2]-d[2]);
    
    std::vector<double> s(3);
    s[0] = (b[0]-d[0]);
    s[1] = (b[1]-d[1]);
    s[2] = (b[2]-d[2]);
    
    std::vector<double> r(3);
    r[0] = (c[0]-d[0]);
    r[1] = (c[1]-d[1]);
    r[2] = (c[2]-d[2]);
    
    
    std::vector<double> cross(3);

        //Cross cross formula
    cross[0] = (s[1] * r[2]) - (s[2] * r[1]);
    cross[1] = (s[2] * r[0]) - (s[0] * r[2]);
    cross[2] = (s[0] * r[1]) - (s[1] * r[0]);
    
    double V = fabs(t[0]*cross[0]+t[1]*cross[1]+t[2]*cross[2])/6.0;
    

    
    return V;
}



double ComputeVolumePrismCell(std::vector<double> P)
{
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    std::vector<double> v2(3);
    std::vector<double> v3(3);
    std::vector<double> v4(3);
    std::vector<double> v5(3);
    
    v0[0] = P[0*3+0]; v1[0] = P[1*3+0];
    v0[1] = P[0*3+1]; v1[1] = P[1*3+1];
    v0[2] = P[0*3+2]; v1[2] = P[1*3+2];
    
    v2[0] = P[2*3+0]; v3[0] = P[3*3+0];
    v2[1] = P[2*3+1]; v3[1] = P[3*3+1];
    v2[2] = P[2*3+2]; v3[2] = P[3*3+2];
    
    v4[0] = P[4*3+0]; v5[0] = P[5*3+0];
    v4[1] = P[4*3+1]; v5[1] = P[5*3+1];
    v4[2] = P[4*3+2]; v5[2] = P[5*3+2];
    
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
    

    
    return V;
    
}





// This function outputs J as an array of 9 values where the matrix is defined as:

/*
 Jac = [J[0], J[1], J[2]
        J[3], J[4], J[5]
        J[6], J[7], J[8]]
*/
// J is computed using the 8-point isoparametric mapping for a hex. The 8-point rule should be sufficient since everything is linear anyways.
std::vector<double> ComputeCenterCoord(std::vector<double> P, int np)
{
 
    std::vector<double> V(3);
    
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
    
    V[0] = 0.0;
    V[1] = 0.0;
    V[2] = 0.0;

    for(int i = 0; i < np; i++)
    {
        N[i] = 1.0/8.0*(1+ref[i*3+0]*eta)*(1+ref[i*3+1]*mu)*(1+ref[i*3+2]*ksi);
                        
        V[0] = V[0]+N[i]*P[i*3+0];
        V[1] = V[1]+N[i]*P[i*3+1];
        V[2] = V[2]+N[i]*P[i*3+2];
    }
    
    delete[] ref;
    delete[] N;
    return V;
}


std::vector<double> ComputeCentroidCoord(std::vector<double> P, int np)
{
    std::vector<double> V(3);
    //V = ComputeCenterCoord(P, np);
    V[0] = 0.0;
    V[1] = 0.0;
    V[2] = 0.0;

    for(int i = 0; i < np; i++)
    {

        V[0] = V[0]+P[i*3+0];
        V[1] = V[1]+P[i*3+1];
        V[2] = V[2]+P[i*3+2];
    }

    V[0] = V[0]/np;
    V[1] = V[1]/np;
    V[2] = V[2]/np;
    
    return V;
}


double* ComputeJAtCenter(std::vector<double> P, int np)
{
    //Vert V;
    double * J = new double[9];
    for(int i=0;i<9;i++)
    {
        J[i] = 0.0;
    }
    int * ref = new int[np*3];
        
    // Allocate the arrays for the mapping function and its derivatives.
    std::vector<double> N(np);
    std::vector<double> dNdeta(np);
    std::vector<double> dNdmu(np);
    std::vector<double> dNdksi(np);
    
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


double GetQualityTetrahedra(std::vector<double> P)
{
    double       h1,h2,h3,h4,h5,h6,det,vol,rap,v1,v2,v3,num;
    double       cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
    double       bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
    

    /* average metric */
    double mm[6];
    mm[0] = 1.0;
    mm[1] = 0.0;
    mm[2] = 0.0;
    mm[3] = 1.0;
    mm[4] = 0.0;
    mm[5] = 1.0;
    
    double* a = new double[3];
    a[0] = P[0*3+0];    a[1] = P[0*3+1];    a[2] = P[0*3+2];
    double* b = new double[3];
    b[0] = P[1*3+0];    b[1] = P[1*3+1];    b[2] = P[1*3+2];
    double* c = new double[3];
    c[0] = P[2*3+0];    c[1] = P[2*3+1];    c[2] = P[2*3+2];
    double* d = new double[3];
    d[0] = P[3*3+0];    d[1] = P[3*3+1];    d[2] = P[3*3+2];
    

    /* volume */
    abx = b[0] - a[0];
    aby = b[1] - a[1];
    abz = b[2] - a[2];

    acx = c[0] - a[0];
    acy = c[1] - a[1];
    acz = c[2] - a[2];

    adx = d[0] - a[0];
    ady = d[1] - a[1];
    adz = d[2] - a[2];

    bcx = c[0] - b[0];
    bcy = c[1] - b[1];
    bcz = c[2] - b[2];

    bdx = d[0] - b[0];
    bdy = d[1] - b[1];
    bdz = d[2] - b[2];

    cdx = d[0] - c[0];
    cdy = d[1] - c[1];
    cdz = d[2] - c[2];
    
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    v1  = acy*adz - acz*ady;
    v2  = acz*adx - acx*adz;
    v3  = acx*ady - acy*adx;
    vol = abx * v1 + aby * v2 + abz * v3;
    if ( vol <= 0. )  return 0.0;

    det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
    if ( det < 1.0e-200 )   {
      return 0.0;
    }
    det = sqrt(det) * vol;

    /* edge lengths */
    h1 = mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz
      + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);
    h2 =  mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz
      + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
    h3 = mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz
      + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);
    h4 =  mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz
      + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);
    h5 =  mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz
      + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
    h6 =  mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz
      + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);

    /* quality */
    rap = h1 + h2 + h3 + h4 + h5 + h6;
    num = sqrt(rap) * rap;

    cal = det / num;


    
    return cal;
}

Array<double>* ComputeJAtCenter_tet_v2(std::vector<double> P)
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
double ComputeDeterminantJ_tet_v2(std::vector<double> P)
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

double ComputeDeterminantJ_tet(std::vector<double> P)
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


double ComputeDeterminantJ(std::vector<double> P, int np)
{
    double* JP1 = ComputeJAtCenter(P, np);
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
    -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
    +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    return DetJ;
    
}



Array<double>* ComputeDeterminantofJacobian(ParArray<int>* ien, Array<double>* xcn)
{
    
    int np = 8; // Assuming its all hexes.
    int nel_loc = ien->getNrow();
    int ncol = ien->getNcol();
    Array<double>* detJ = new Array<double>(nel_loc,1);
    //double* P = new double[np*3];
    std::vector<double> P(np*3);
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
    //delete[] P;
    
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
    //double* P = new double[np*3];
    std::vector<double> P(np*3);
    
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
    
    //delete[] P;
    
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
    //double* P = new double[np*3];
    std::vector<double> P(np*3);
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
    //delete[] P;
    
    return vol_verts;
}



void UnitTestJacobian()
{
    //double* Hex = new double[8*3];
    std::vector<double> Hex(8*3);
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



std::map<int,std::vector<std::vector<double> > > ComputeMetric_Lite(MPI_Comm comm, 
                        RepartitionObject* tetra_repart,
                        std::map<int,std::vector<double> > tetra_grad, 
                        std::map<int,std::vector<double> > &eigvalues, 
                        Inputs* inputs)
{
    // Preparing the metric tensor field.
    std::vector<double> eignval(3,0.0);
    double po = 6.0;

    int rec         = inputs->recursive;
    int ext         = inputs->extended;
    double hmin     = inputs->hmin;
    double hmax     = inputs->hmax;
    double Scale    = inputs->MetScale;
    std::map<int,std::vector<std::vector<double> > > metric_vmap;



    std::map<int,std::vector<double> >::iterator itmidv;
    int yep = 0;
    std::vector<int> Owned_Elem_t                           = tetra_repart->getLocElem();
    std::map<int,std::vector<int> > gE2gV_t                 = tetra_repart->getElement2VertexMap();
    std::map<int,std::set<int> > node2element_map    = tetra_repart->GetNode2ElementMap();

    for(int i=0;i<Owned_Elem_t.size();i++)
    {
        int elid = Owned_Elem_t[i];
        int nv   = gE2gV_t[elid].size();
        for(int j=0;j<nv;j++)
        {
            int gvid = gE2gV_t[elid][j];

            if(node2element_map.find(gvid)!=node2element_map.end())
            {
                std::set<int>::iterator its;
                std::set<int> elems = node2element_map[gvid];
                double gval         = 0.0;
                int nc              = 0;
                std::vector<double> row_tmp(6,0.0);

                for(its=elems.begin();its!=elems.end();its++)
                {
                    int elid = *its;

                    if(tetra_grad.find(elid)!=tetra_grad.end())
                    {
                        for(int k=0;k<6;k++)
                        {
                            row_tmp[k]=row_tmp[k]+tetra_grad[elid][3+k];
                        }
                        nc++;
                    }   
                }
                
                for(int k=0;k<6;k++)
                {
                    row_tmp[k]  = row_tmp[k]/nc;
                }

                std::vector<double> row(9,0.0);

                row[0]=row_tmp[0];  row[1]=row_tmp[1];  row[2]=row_tmp[2];
                row[3]=row_tmp[1];  row[4]=row_tmp[3];  row[5]=row_tmp[4];
                row[6]=row_tmp[2];  row[7]=row_tmp[4];  row[8]=row_tmp[5];
                
                Eig* eig = ComputeEigenDecomp(3, row.data());
                std::vector<double> Dre(3,0.0);
                Dre[0] = eig->Dre[0];
                Dre[1] = eig->Dre[1];
                Dre[2] = eig->Dre[2];
                eigvalues[elid] = Dre;

                eignval[0] =  std::min(std::max(Scale*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
                eignval[1] =  std::min(std::max(Scale*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
                eignval[2] =  std::min(std::max(Scale*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
                
                double Lambdamax = *std::max_element(eignval.begin(),eignval.end());
                double Lambdamin = *std::min_element(eignval.begin(),eignval.end());
                //std::cout << "eign " << eignval[0] << " " << eignval[1] << " " << eignval[2] << " --> " << eig->Dre[0] << " " << eig->Dre[1] << " " << eig->Dre[2] << std::endl;
                std::vector<std::vector<double> > Diag(3);
                std::vector<std::vector<double> > EigVec(3);
                
                for(int k=0;k<3;k++)
                {
                    std::vector<double> rowDiag(3,0.0);
                    rowDiag[k]      = eignval[k];
                    std::vector<double> rowEigVec(3,0.0);
                    Diag[k]         = rowDiag;
                    EigVec[k]       = rowEigVec;
                }

                EigVec[0][0] = eig->V[0];   EigVec[0][1] = eig->V[1];   EigVec[0][2] = eig->V[2];
                EigVec[1][0] = eig->V[3];   EigVec[1][1] = eig->V[4];   EigVec[1][2] = eig->V[5];
                EigVec[2][0] = eig->V[6];   EigVec[2][1] = eig->V[7];   EigVec[2][2] = eig->V[8];

                std::vector<std::vector<double> > iVR       = MatInv_Lite(EigVec);
                std::vector<std::vector<double> > Rs        = MatMul_Lite(EigVec,Diag);  
                std::vector<std::vector<double> > metric    = MatMul_Lite(Rs,iVR);

                double detMetric        = metric[0][0]*(metric[1][1]*metric[2][2]-metric[2][1]*metric[1][2])
                                        - metric[0][1]*(metric[1][0]*metric[2][2]-metric[2][0]*metric[1][2])
                                        + metric[0][2]*(metric[1][0]*metric[2][1]-metric[2][0]*metric[1][1]);
                
                double pow              = -1.0/(2.0*po+3.0);
                double eigRat           = Lambdamin/Lambdamax;
                detMetric               = std::pow(detMetric,pow);

                for(int i=0;i<3;i++)
                {
                    for(int j=0;j<3;j++)
                    {
                        metric[i][j] = detMetric*metric[i][j];
                    }
                }

                metric_vmap[gvid] = metric;

                delete[] eig->Dre;
                delete[] eig->Dim;
                delete[] eig->V;
                delete[] eig->iV;
            }
        }
    }



    //==================================================================================================
    // std::map<int,std::vector<double> > metric_diagnose;
    // for(int i=0;i<Owned_Elem_t.size();i++)
    // {
    //     int elid = Owned_Elem_t[i];

    //     std::vector<double> metric_tmp(9,0.0);
    //     std::vector<double> row_tmp(6,0.0);
    //     for(int k=0;k<6;k++)
    //     {
    //         row_tmp[k]=row_tmp[k]+tetra_grad[elid][3+k];
    //     }


    //     metric_tmp[0]=row_tmp[0];  metric_tmp[1]=row_tmp[1];  metric_tmp[2]=row_tmp[2];
    //     metric_tmp[3]=row_tmp[1];  metric_tmp[4]=row_tmp[3];  metric_tmp[5]=row_tmp[4];
    //     metric_tmp[6]=row_tmp[2];  metric_tmp[7]=row_tmp[4];  metric_tmp[8]=row_tmp[5];

    
    //     Eig* eig = ComputeEigenDecomp(3, metric_tmp.data());
    //     eignval[0] =  std::min(std::max(fabs(Scale*eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
    //     eignval[1] =  std::min(std::max(fabs(Scale*eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
    //     eignval[2] =  std::min(std::max(fabs(Scale*eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
    //     double Lambdamax = *std::max_element(eignval.begin(),eignval.end());
    //     double Lambdamin = *std::min_element(eignval.begin(),eignval.end());

    //     std::vector<std::vector<double> > Diag(3);
    //     std::vector<std::vector<double> > EigVec(3);
        
    //     for(int k=0;k<3;k++)
    //     {
    //         std::vector<double> rowDiag(3,0.0);
    //         rowDiag[k] = eignval[k];
    //         std::vector<double> rowEigVec(3,0.0);
    //         Diag[k]      = rowDiag;
    //         EigVec[k]    = rowEigVec;
    //     }
    
    //     EigVec[0][0] = eig->V[0];   EigVec[0][1] = eig->V[1];   EigVec[0][2] = eig->V[2];
    //     EigVec[1][0] = eig->V[3];   EigVec[1][1] = eig->V[4];   EigVec[1][2] = eig->V[5];
    //     EigVec[2][0] = eig->V[6];   EigVec[2][1] = eig->V[7];   EigVec[2][2] = eig->V[8];

    //     std::vector<std::vector<double> > iVR       = MatInv_Lite(EigVec);
    //     std::vector<std::vector<double> > Rs        = MatMul_Lite(EigVec,Diag);  
    //     std::vector<std::vector<double> > metric    = MatMul_Lite(Rs,iVR);

    //     double detMetric        = metric[0][0]*(metric[1][1]*metric[2][2]-metric[2][1]*metric[1][2])
    //                             - metric[0][1]*(metric[1][0]*metric[2][2]-metric[2][0]*metric[1][2])
    //                             + metric[0][2]*(metric[1][0]*metric[2][1]-metric[2][0]*metric[1][1]);
        
    //     double pow              = -1.0/(2.0*po+3.0);
    //     double eigRat           = Lambdamin/Lambdamax;
    //     detMetric               = std::pow(detMetric,pow);

    //     for(int i=0;i<3;i++)
    //     {
    //         for(int j=0;j<3;j++)
    //         {
    //             metric_tmp[i*3+j] = detMetric*metric[i][j];
    //         }
    //     }

    //     metric_diagnose[elid] = metric_tmp;
    // }

    // string filename_metric = "metric_" + std::to_string(world_rank) + ".vtu";
    // std::map<int,std::string > varnames_metric;
    // varnames_metric[0] = "m00";    varnames_metric[1] = "m01";    varnames_metric[2] = "m02";
    // varnames_metric[3] = "m10";    varnames_metric[4] = "m11";    varnames_metric[5] = "m12";
    // varnames_metric[6] = "m20";    varnames_metric[7] = "m21";    varnames_metric[8] = "m22";

    // OutputTetraMeshPartitionVTK(comm,
    //                         filename_metric, 
    //                         Owned_Elem_t, 
    //                         gE2gV_t, 
    //                         metric_diagnose, 
    //                         varnames_metric, 
    //                         LocalVertsMap_t);
    //==================================================================================================

    return metric_vmap;

}





void ComputeMetricWithWake(Partition* Pa,
                   MPI_Comm comm,
                   std::map<int,Array<double>* > scale_vm,
                   std::map<int,Array<double>* > &Hess_vm,
                   double sumvol, double po, double hwake, int recursive, double hmin, double hmax, double MetScale)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++

    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double* Hmet = new double[9];

    std::map<int,Array<double>*>::iterator itm;
    int i = 0;
    double Lambdamax,Lambdamin;
    std::vector<double> eignval(3);
    int anitel = 0;
    for(itm=Hess_vm.begin();itm!=Hess_vm.end();itm++)
    {
        int glob_vid = itm->first;

        if(recursive == 0)
        {
            Hmet[0] = Hess_vm[glob_vid]->getVal(3,0);
            Hmet[1] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[2] = Hess_vm[glob_vid]->getVal(5,0);

            Hmet[3] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[4] = Hess_vm[glob_vid]->getVal(6,0);
            Hmet[5] = Hess_vm[glob_vid]->getVal(7,0);

            Hmet[6] = Hess_vm[glob_vid]->getVal(5,0);
            Hmet[7] = Hess_vm[glob_vid]->getVal(7,0);
            Hmet[8] = Hess_vm[glob_vid]->getVal(8,0);
        }
        if(recursive == 1)
        {
            Hmet[0] = Hess_vm[glob_vid]->getVal(0,0);
            Hmet[1] = Hess_vm[glob_vid]->getVal(1,0);
            Hmet[2] = Hess_vm[glob_vid]->getVal(2,0);

            Hmet[3] = Hess_vm[glob_vid]->getVal(1,0);
            Hmet[4] = Hess_vm[glob_vid]->getVal(3,0);
            Hmet[5] = Hess_vm[glob_vid]->getVal(4,0);

            Hmet[6] = Hess_vm[glob_vid]->getVal(2,0);
            Hmet[7] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[8] = Hess_vm[glob_vid]->getVal(5,0);
        }
        
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
        
        eignval[0] =  std::min(std::max(MetScale*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[1] =  std::min(std::max(MetScale*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[2] =  std::min(std::max(MetScale*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        Lambdamax = *std::max_element(eignval.begin(),eignval.end());
        Lambdamin = *std::min_element(eignval.begin(),eignval.end());
        
        DR->setVal(0,0,eignval[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,eignval[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,eignval[2]);
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        
        double detRf = Rf->getVal(0,0)*(Rf->getVal(1,1)*Rf->getVal(2,2)-Rf->getVal(2,1)*Rf->getVal(1,2))
                      -Rf->getVal(0,1)*(Rf->getVal(1,0)*Rf->getVal(2,2)-Rf->getVal(2,0)*Rf->getVal(1,2))
                      +Rf->getVal(0,2)*(Rf->getVal(1,0)*Rf->getVal(2,1)-Rf->getVal(2,0)*Rf->getVal(1,1));
        
        
        double pow = -1.0/(2.0*po+3.0);
        double eigRat = Lambdamin/Lambdamax;
        
        detRf = std::pow(detRf,pow);
        Array<double>* Habs  = new Array<double>(3,3);
	
        if(scale_vm[glob_vid]->getVal(0,0)>0.1)
        {
            double wi    = scale_vm[glob_vid]->getVal(0,0);
            double wa    = 1.0-scale_vm[glob_vid]->getVal(0,0);
//            
//            wi = 0.0;
//            wa = 1.0;
            double hiso  = hwake;
           
            //std::cout << wi << " " << wa << std::endl;
	    
            Habs->setVal(0,0,wi/(hiso*hiso)+wa*sumvol*detRf*Rf->getVal(0,0));
            Habs->setVal(0,1,wa*sumvol*detRf*Rf->getVal(0,1));
            Habs->setVal(0,2,wa*sumvol*detRf*Rf->getVal(0,2));
            
            Habs->setVal(1,0,wa*sumvol*detRf*Rf->getVal(1,0));
            Habs->setVal(1,1,wi/(hiso*hiso)+wa*sumvol*detRf*Rf->getVal(1,1));
            Habs->setVal(1,2,wa*sumvol*detRf*Rf->getVal(1,2));

            Habs->setVal(2,0,wa*sumvol*detRf*Rf->getVal(2,0));
            Habs->setVal(2,1,wa*sumvol*detRf*Rf->getVal(2,1));
            Habs->setVal(2,2,wi/(hiso*hiso)+wa*sumvol*detRf*Rf->getVal(2,2));
            
            anitel++;
        }
        else
        {
            Habs->setVal(0,0,sumvol*detRf*Rf->getVal(0,0));
            Habs->setVal(0,1,sumvol*detRf*Rf->getVal(0,1));
            Habs->setVal(0,2,sumvol*detRf*Rf->getVal(0,2));
            
            Habs->setVal(1,0,sumvol*detRf*Rf->getVal(1,0));
            Habs->setVal(1,1,sumvol*detRf*Rf->getVal(1,1));
            Habs->setVal(1,2,sumvol*detRf*Rf->getVal(1,2));

            Habs->setVal(2,0,sumvol*detRf*Rf->getVal(2,0));
            Habs->setVal(2,1,sumvol*detRf*Rf->getVal(2,1));
            Habs->setVal(2,2,sumvol*detRf*Rf->getVal(2,2));
        }
        
        Hess_vm[glob_vid]=Habs;
        
        delete Rf;
        delete DR;
        delete UR;
        delete Rs;
        //delete Rf;

        delete[] iV;
        delete[] V;
        delete[] WR;
        delete[] WI;
        i++;
    }
    
    int anitel_red;
    MPI_Allreduce(&anitel, &anitel_red, 1, MPI_INT, MPI_SUM, comm);

    delete[] Hmet;
}










std::map<int,std::vector<std::vector<double> > > ComputeElementMetric_Lite(MPI_Comm comm, 
                        RepartitionObject* tetra_repart,
                        std::map<int,std::vector<double> > tetra_grad,
                        std::map<int,std::vector<double> > &eigvalues, 
                        Inputs* inputs)
{
    // Preparing the metric tensor field.
    std::vector<double> eignval(3,0.0);
    double po = 6.0;

    int rec         = inputs->recursive;
    int ext         = inputs->extended;
    double hmin     = inputs->hmin;
    double hmax     = inputs->hmax;
    double Scale    = inputs->MetScale;
    std::map<int,std::vector<std::vector<double> > > metric_vmap;

    std::map<int,std::vector<double> >::iterator itmidv;
    int yep = 0;
    std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();
    std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
    std::map<int,std::set<int> > node2element_map       = tetra_repart->GetNode2ElementMap();

    for(int i=0;i<Owned_Elem_t.size();i++)
    {
        int elid = Owned_Elem_t[i];

        std::vector<double> row(9,0.0);

        row[0]=tetra_grad[elid][3];  row[1]=tetra_grad[elid][4];  row[2]=tetra_grad[elid][5];
        row[3]=tetra_grad[elid][5];  row[4]=tetra_grad[elid][6];  row[5]=tetra_grad[elid][7];
        row[6]=tetra_grad[elid][5];  row[7]=tetra_grad[elid][7];  row[8]=tetra_grad[elid][8];
        
        Eig* eig = ComputeEigenDecomp(3, row.data());

        std::vector<double> Dre(3,0.0);
        Dre[0] = eig->Dre[0];
        Dre[1] = eig->Dre[1];
        Dre[2] = eig->Dre[2];
        eigvalues[elid] = Dre;

        eignval[0] =  std::min(std::max(Scale*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[1] =  std::min(std::max(Scale*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[2] =  std::min(std::max(Scale*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        double Lambdamax = *std::max_element(eignval.begin(),eignval.end());
        double Lambdamin = *std::min_element(eignval.begin(),eignval.end());
        //std::cout << "eign " << eignval[0] << " " << eignval[1] << " " << eignval[2] << " --> " << eig->Dre[0] << " " << eig->Dre[1] << " " << eig->Dre[2] << std::endl;
        std::vector<std::vector<double> > Diag(3);
        std::vector<std::vector<double> > EigVec(3);
        
        for(int k=0;k<3;k++)
        {
            std::vector<double> rowDiag(3,0.0);
            rowDiag[k]      = eignval[k];
            std::vector<double> rowEigVec(3,0.0);
            Diag[k]         = rowDiag;
            EigVec[k]       = rowEigVec;
        }
    
        EigVec[0][0] = eig->V[0];   EigVec[0][1] = eig->V[1];   EigVec[0][2] = eig->V[2];
        EigVec[1][0] = eig->V[3];   EigVec[1][1] = eig->V[4];   EigVec[1][2] = eig->V[5];
        EigVec[2][0] = eig->V[6];   EigVec[2][1] = eig->V[7];   EigVec[2][2] = eig->V[8];

        std::vector<std::vector<double> > iVR       = MatInv_Lite(EigVec);
        std::vector<std::vector<double> > Rs        = MatMul_Lite(EigVec,Diag);  
        std::vector<std::vector<double> > metric    = MatMul_Lite(Rs,iVR);

        double detMetric        = metric[0][0]*(metric[1][1]*metric[2][2]-metric[2][1]*metric[1][2])
                                - metric[0][1]*(metric[1][0]*metric[2][2]-metric[2][0]*metric[1][2])
                                + metric[0][2]*(metric[1][0]*metric[2][1]-metric[2][0]*metric[1][1]);
        
        double pow              = -1.0/(2.0*po+3.0);
        double eigRat           = Lambdamin/Lambdamax;
        detMetric               = std::pow(detMetric,pow);

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                metric[i][j] = detMetric*metric[i][j];
            }
        }

        metric_vmap[elid] = metric;

        delete[] eig->Dre;
        delete[] eig->Dim;
        delete[] eig->V;
        delete[] eig->iV;
    }

    return metric_vmap;
}









void ComputeMetric(Partition* Pa, MPI_Comm comm, std::map<int,Array<double>* > &Hess_vm,
                   double sumvol, double po, int recursive, int extended, double hmin, double hmax, double MetScale)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++
    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double* Hmet = new double[9];

    std::map<int,Array<double>*>::iterator itm;
    int i = 0;
    double fiso = 0.01;
    double Lambdamax,Lambdamin;
    std::vector<double> eignval(3);
    int anitel = 0;
    for(itm=Hess_vm.begin();itm!=Hess_vm.end();itm++)
    {
        int glob_vid = itm->first;
        
        if(Hess_vm[glob_vid]->getNrow()==9)
        {
	    //std::cout << "Hexx" << glob_vid << " " << Hess_vm[glob_vid]->getVal(3,0) << " "  << Hess_vm[glob_vid]->getVal(3,0) << " " << Hess_vm[glob_vid]->getVal(5,0) << " " << Hess_vm[glob_vid]->getVal(7,0) << " " << Hess_vm[glob_vid]->getVal(8,0) << std::endl; 
            Hmet[0] = Hess_vm[glob_vid]->getVal(3,0);
            Hmet[1] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[2] = Hess_vm[glob_vid]->getVal(5,0);

            Hmet[3] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[4] = Hess_vm[glob_vid]->getVal(6,0);
            Hmet[5] = Hess_vm[glob_vid]->getVal(7,0);

            Hmet[6] = Hess_vm[glob_vid]->getVal(5,0);
            Hmet[7] = Hess_vm[glob_vid]->getVal(7,0);
            Hmet[8] = Hess_vm[glob_vid]->getVal(8,0);
        }
        if(Hess_vm[glob_vid]->getNrow()==6)
        {
            Hmet[0] = Hess_vm[glob_vid]->getVal(0,0);
            Hmet[1] = Hess_vm[glob_vid]->getVal(1,0);
            Hmet[2] = Hess_vm[glob_vid]->getVal(2,0);

            Hmet[3] = Hess_vm[glob_vid]->getVal(1,0);
            Hmet[4] = Hess_vm[glob_vid]->getVal(3,0);
            Hmet[5] = Hess_vm[glob_vid]->getVal(4,0);

            Hmet[6] = Hess_vm[glob_vid]->getVal(2,0);
            Hmet[7] = Hess_vm[glob_vid]->getVal(4,0);
            Hmet[8] = Hess_vm[glob_vid]->getVal(5,0);
        }
        
        
        
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
        
        eignval[0] =  std::min(std::max(MetScale*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[1] =  std::min(std::max(MetScale*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        eignval[2] =  std::min(std::max(MetScale*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        Lambdamax = *std::max_element(eignval.begin(),eignval.end());
        Lambdamin = *std::min_element(eignval.begin(),eignval.end());
        
        DR->setVal(0,0,eignval[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,eignval[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,eignval[2]);
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        
        double detRf = Rf->getVal(0,0)*(Rf->getVal(1,1)*Rf->getVal(2,2)-Rf->getVal(2,1)*Rf->getVal(1,2))
                      -Rf->getVal(0,1)*(Rf->getVal(1,0)*Rf->getVal(2,2)-Rf->getVal(2,0)*Rf->getVal(1,2))
                      +Rf->getVal(0,2)*(Rf->getVal(1,0)*Rf->getVal(2,1)-Rf->getVal(2,0)*Rf->getVal(1,1));
        
        
        double pow = -1.0/(2.0*po+3.0);
        double eigRat = Lambdamin/Lambdamax;
        
        detRf = std::pow(detRf,pow);
        Array<double>* Habs  = new Array<double>(3,3);

        Habs->setVal(0,0,sumvol*detRf*Rf->getVal(0,0));
        Habs->setVal(0,1,sumvol*detRf*Rf->getVal(0,1));
        Habs->setVal(0,2,sumvol*detRf*Rf->getVal(0,2));
        
        Habs->setVal(1,0,sumvol*detRf*Rf->getVal(1,0));
        Habs->setVal(1,1,sumvol*detRf*Rf->getVal(1,1));
        Habs->setVal(1,2,sumvol*detRf*Rf->getVal(1,2));

        Habs->setVal(2,0,sumvol*detRf*Rf->getVal(2,0));
        Habs->setVal(2,1,sumvol*detRf*Rf->getVal(2,1));
        Habs->setVal(2,2,sumvol*detRf*Rf->getVal(2,2));
/*     
        Array<double>* Habs  = new Array<double>(3,3);
	Habs->setVal(0,0,f*1.0);
        Habs->setVal(0,1,0.0);
        Habs->setVal(0,2,0.0);
        
        Habs->setVal(1,0,0.0);
        Habs->setVal(1,1,f*1.0);
        Habs->setVal(1,2,0.0);

        Habs->setVal(2,0,0.0);
        Habs->setVal(2,1,0.0);
        Habs->setVal(2,2,f*1.0);
*/        
        //std::cout << "Habs " <<  Habs->getVal(0,0) << " " << Habs->getVal(0,1) << " " << Habs->getVal(0,2) << " " << Habs->getVal(1,1) << " " << Habs->getVal(1,2) << " " << Habs->getVal(2,2) << std::endl; 
        Hess_vm[glob_vid]=Habs;
        
   //     delete Rf;
        delete DR;
        delete UR;
    //    delete Rs;
        //delete Rf;

        delete[] iV;
        delete[] V;
        delete[] WR;
        delete[] WI;
        
        i++;
    }
    int anitel_red;
    MPI_Allreduce(&anitel, &anitel_red, 1, MPI_INT, MPI_SUM, comm);
    
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
    
    std::map<int,Array<double>* > Uelem_all;// = P->PartitionAuxilaryData(U, comm);
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
            u_o      = Uelem_all[adjEl_id]->getVal(0,0);
            
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
    std::map<int,int> gE2lE                    = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem                  = Pa->getLocElem();
    int nLocElem                               = Loc_Elem.size();
    std::map<int,std::vector<int> > gE2lV      = Pa->getGlobElem2LocVerts();
    std::vector<std::vector<double> > locVerts = Pa->getLocalVerts();
    Array<double>* Volumes = new Array<double>(nLocElem,1);
    std::vector<int> vijkIDs;
    //double* Pijk = new double[8*3];
    std::vector<double> Pijk(8*3);
    for(i=0;i<nLocElem;i++)
    {
       int gEl = Loc_Elem[i];

       vijkIDs = gE2lV[gEl];

       for(k=0;k<vijkIDs.size();k++)
       {
          loc_vid     = vijkIDs[k];
          Pijk[k*3+0] = locVerts[loc_vid][0];
          Pijk[k*3+1] = locVerts[loc_vid][1];
          Pijk[k*3+2] = locVerts[loc_vid][2];
       }

       double Vol = ComputeVolumeHexCell(Pijk);
       Volumes->setVal(i,0,Vol);
    }
    
    //delete[] Pijk;
    return Volumes;
}


