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


inline double ComputeEdgeLength(Vert* v0, Vert* v1)
{
    return sqrt((v0->x - v1->x) * (v0->x - v1->x)+
                (v0->y - v1->y) * (v0->y - v1->y)+
                (v0->z - v1->z) * (v0->z - v1->z));
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
    
    Vert v0;
    Vert v1;
    
    v0.x = P[0*3+0]; v1.x = P[1*3+0];
    v0.y = P[0*3+1]; v1.y = P[1*3+1];
    v0.z = P[0*3+2]; v1.z = P[1*3+2];
    
    L01 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[1*3+0]; v1.x = P[5*3+0];
    v0.y = P[1*3+1]; v1.y = P[5*3+1];
    v0.z = P[1*3+2]; v1.z = P[5*3+2];
    
    L15 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[1*3+0]; v1.x = P[2*3+0];
    v0.y = P[1*3+1]; v1.y = P[2*3+1];
    v0.z = P[1*3+2]; v1.z = P[2*3+2];
    
    H12 = ComputeEdgeLength(&v0,&v1);
    b0 = 0.5*L01*L15;
    double vol0 = 1.0/3.0*b0*H12;
    //==================================================

    v0.x = P[0*3+0]; v1.x = P[4*3+0];
    v0.y = P[0*3+1]; v1.y = P[4*3+1];
    v0.z = P[0*3+2]; v1.z = P[4*3+2];
    
    L04 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[4*3+0]; v1.x = P[5*3+0];
    v0.y = P[4*3+1]; v1.y = P[5*3+1];
    v0.z = P[4*3+2]; v1.z = P[5*3+2];
    
    L45 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[4*3+0]; v1.x = P[7*3+0];
    v0.y = P[4*3+1]; v1.y = P[7*3+1];
    v0.z = P[4*3+2]; v1.z = P[7*3+2];
    
    H47 = ComputeEdgeLength(&v0,&v1);
    b1 = 0.5*L04*L45;
    double vol1 = 1.0/3.0*b1*H47;
    
    //==================================================
    
    v0.x = P[3*3+0]; v1.x = P[7*3+0];
    v0.y = P[3*3+1]; v1.y = P[7*3+1];
    v0.z = P[3*3+2]; v1.z = P[7*3+2];
    
    L37 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[2*3+0]; v1.x = P[3*3+0];
    v0.y = P[2*3+1]; v1.y = P[3*3+1];
    v0.z = P[2*3+2]; v1.z = P[3*3+2];
    
    L23 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[3*3+0]; v1.x = P[0*3+0];
    v0.y = P[3*3+1]; v1.y = P[0*3+1];
    v0.z = P[3*3+2]; v1.z = P[0*3+2];
   
    H30 = ComputeEdgeLength(&v0,&v1);
    b2 = 0.5*L37*L23;
    double vol2 = 1.0/3.0*b2*H30;
    
    //==================================================
    
    v0.x = P[2*3+0]; v1.x = P[6*3+0];
    v0.y = P[2*3+1]; v1.y = P[6*3+1];
    v0.z = P[2*3+2]; v1.z = P[6*3+2];
    
    L26 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[6*3+0]; v1.x = P[7*3+0];
    v0.y = P[6*3+1]; v1.y = P[7*3+1];
    v0.z = P[6*3+2]; v1.z = P[7*3+2];
    
    L67 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[5*3+0]; v1.x = P[6*3+0];
    v0.y = P[5*3+1]; v1.y = P[6*3+1];
    v0.z = P[5*3+2]; v1.z = P[6*3+2];
    
    H56 = ComputeEdgeLength(&v0,&v1);
    b3 = 0.5*L26*L67;
    double vol3 = 1.0/3.0*b3*H56;
    
    return vol0+vol1+vol2+vol3;
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

Array<double>* ComputeMetric(std::vector<Vert> Verts, Array<double>* grad, Array<double>* hessian, double max_v, std::vector<std::vector<int> > loc_elem2verts_loc, int nloc, MPI_Comm comm, int dim)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++Required Parameter set for scaling eigenvalues/eigenvectors+++++++++++++++
    double R            = 0.056;
    double hmin         = 0.00001;
    double hmax         = 0.1;
    
    double f            = 0.15;
    std::cout << "Metric Tensor Field gets computed..." << std::endl;
    //double d2udx2_v_max = *std::max_element(d2udx2_v.begin(), d2udx2_v.end());
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double* Hmet;
    if(dim==3)
    {
        Hmet = new double[9];
    }
    if(dim==2)
    {
        Hmet = new double[4];
    }
    int nVerts = Verts.size();
    int i;
    Array<double>* metric = new Array<double>(nVerts,9);
//    string filename11 = "metric_plane.dat";
//    ofstream myfile11;
    std::string filename = "metricTensor_rank_" + std::to_string(rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\",  \"drhodx\",  \"drhody\",  \"drhodz\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\", \"H00\", \"H01\", \"H02\", \"H11\", \"H12\", \"H22\"" << std::endl;
    myfile <<"ZONE N = " << nVerts << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    for(i=0;i<nVerts;i++)
    {
        if(dim == 2)
        {
            Hmet[0] = hessian->getVal(i,0);
            Hmet[1] = hessian->getVal(i,1);
            Hmet[2] = hessian->getVal(i,1);
            Hmet[3] = hessian->getVal(i,3);
        
            double * WR = new double[3];
            double * WI = new double[3];
            double * V  = new double[3*3];
            double * iV = new double[3*3];
            double* WRn = new double[3];
            
            Array<double>* DR  = new Array<double>(3,3);
            Array<double>* VR  = new Array<double>(3,3);
            Array<double>* UR  = new Array<double>(3,3);
            for(int j=0;j<3;j++)
            {
                WR[j]  = 0.0;
                WRn[j] = 0.0;
                WI[j]  = 0.0;
                
                for(int k=0;k<3;k++)
                {
                    DR->setVal(j,k,0.0);
                    VR->setVal(j,k,0.0);
                    //iVR->setVal(j,k,0.0);
                    V[j*3+k] = 0.0;
                    iV[j*3+k] = 0.0;
                }
            }
            SVD* svd = ComputeSVD(2,2,Hmet);
            Eig* eig = ComputeEigenDecomp(2, Hmet);
            
            WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            WRn[2] = 1.0/(hmax*hmax);
            
            DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
            DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
            DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);

            UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,0.0);
            UR->setVal(1,0,eig->V[2]);UR->setVal(1,1,eig->V[3]);UR->setVal(1,2,0.0);
            UR->setVal(2,0,0.0);      UR->setVal(2,1,0.0);      UR->setVal(2,2,1.0);
            
            Array<double>* iVR = MatInv(UR);
            Array<double>* Rs = MatMul(UR,DR);
            Array<double>* Rf = MatMul(Rs,iVR);
            
            metric->setVal(i,0, Rf->getVal(0,0));
            metric->setVal(i,1, Rf->getVal(0,1));
            metric->setVal(i,2, Rf->getVal(0,2));
            
            metric->setVal(i,3, Rf->getVal(1,0));
            metric->setVal(i,4, Rf->getVal(1,1));
            metric->setVal(i,5, Rf->getVal(1,2));
            
            metric->setVal(i,6, Rf->getVal(2,0));
            metric->setVal(i,7, Rf->getVal(2,1));
            metric->setVal(i,8, Rf->getVal(2,2));
            
            myfile << Verts[i].x << " " << Verts[i].y << " " << Verts[i].z << " " << Rf->getVal(0,0) << " " << Rf->getVal(0,1) << " " << Rf->getVal(0,2) << " " << Rf->getVal(1,1) << " " << Rf->getVal(1,2) << " " << Rf->getVal(2,2) << std::endl;
            
            delete DR;
            delete UR;
            delete VR;
            delete Rs;
            delete Rf;

            delete[] iV;
            delete[] V;
            delete[] WR;
            delete[] WI;
            delete[] WRn;
        }
        
        if(dim == 3)
        {
            Hmet[0] = hessian->getVal(i,0);
            Hmet[1] = hessian->getVal(i,1);
            Hmet[2] = hessian->getVal(i,2);

            Hmet[3] = hessian->getVal(i,1);
            Hmet[4] = hessian->getVal(i,3);
            Hmet[5] = hessian->getVal(i,4);

            Hmet[6] = hessian->getVal(i,2);
            Hmet[7] = hessian->getVal(i,4);
            Hmet[8] = hessian->getVal(i,5);
/*
	    Hmet[0] = hessian->getVal(i,0);
            Hmet[1] = hessian->getVal(i,1);
            Hmet[2] = 1.0e-06;

            Hmet[3] = hessian->getVal(i,1);
            Hmet[4] = hessian->getVal(i,3);
            Hmet[5] = 1.0e-06;

            Hmet[6] = 1.0e-06;
            Hmet[7] = 1.0e-06;
            Hmet[8] = 1.0e-06;
*/
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
            SVD* svd = ComputeSVD(3,3,Hmet);
            Eig* eig = ComputeEigenDecomp(3, Hmet);
            
//            WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
//            WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
//            WRn[2] = std::min(std::max(f*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            
            WRn[0] = std::min(std::max(f*fabs(svd->s[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            WRn[1] = std::min(std::max(f*fabs(svd->s[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            WRn[2] = std::min(std::max(f*fabs(svd->s[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            
            DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
            DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
            DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);

//            UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
//            UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
//            UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);
            
//            UR->setVal(0,0,svd->u[0]);UR->setVal(0,1,svd->u[1]);UR->setVal(0,2,svd->u[2]);
//            UR->setVal(1,0,svd->u[3]);UR->setVal(1,1,svd->u[4]);UR->setVal(1,2,svd->u[5]);
//            UR->setVal(2,0,svd->u[6]);UR->setVal(2,1,svd->u[7]);UR->setVal(2,2,svd->u[8]);
            
            UR->setVal(0,0,svd->vt[0]);UR->setVal(0,1,svd->vt[1]);UR->setVal(0,2,svd->vt[2]);
            UR->setVal(1,0,svd->vt[3]);UR->setVal(1,1,svd->vt[4]);UR->setVal(1,2,svd->vt[5]);
            UR->setVal(2,0,svd->vt[6]);UR->setVal(2,1,svd->vt[7]);UR->setVal(2,2,svd->vt[8]);
            
            Array<double>* iVR = MatInv(UR);
            Array<double>* Rs = MatMul(UR,DR);
            Array<double>* Rf = MatMul(Rs,iVR);
            
            metric->setVal(i,0, Rf->getVal(0,0));
            metric->setVal(i,1, Rf->getVal(0,1));
            metric->setVal(i,2, Rf->getVal(0,2));
            
            metric->setVal(i,3, Rf->getVal(1,0));
            metric->setVal(i,4, Rf->getVal(1,1));
            metric->setVal(i,5, Rf->getVal(1,2));
            
            metric->setVal(i,6, Rf->getVal(2,0));
            metric->setVal(i,7, Rf->getVal(2,1));
            metric->setVal(i,8, Rf->getVal(2,2));
            
            myfile << Verts[i].x << " " << Verts[i].y << " " << Verts[i].z << " " << Rf->getVal(0,0) << " " << Rf->getVal(0,1) << " " << Rf->getVal(0,2) << " " << Rf->getVal(1,1) << " " << Rf->getVal(1,2) << " " << Rf->getVal(2,2) << " " << hessian->getVal(i,0)<< " " << hessian->getVal(i,1)<< " " << hessian->getVal(i,2)<< " " << hessian->getVal(i,3)<< " " << hessian->getVal(i,4)<< " " << hessian->getVal(i,5) << std::endl;
            
            delete DR;
            delete UR;
            delete Rs;
            delete Rf;

            delete[] iV;
            delete[] V;
            delete[] WR;
            delete[] WI;
            delete[] WRn;
        }
        
    }
    
    
    for(int i=0;i<nloc;i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    
    myfile.close();
    
    
    
    delete[] Hmet;
    return metric;
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
    std::vector<Vert> locVerts              = Pa->getLocalVerts();
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
          Pijk[k*3+0] = locVerts[loc_vid].x;
          Pijk[k*3+1] = locVerts[loc_vid].y;
          Pijk[k*3+2] = locVerts[loc_vid].z;
       }

       double Vol = ComputeVolumeHexCell(Pijk);
       Volumes->setVal(i,0,Vol);
    }
    
    delete[] Pijk;
    return Volumes;
}



/*
double* ComputeVolumeCells(Array<double>* xcn, Array<int>* ien)
{
    int Nelements = ien->nrow;
    double * vol_cells = new double[Nelements];
    int np = 8;
    double* P = new double[np*3];
    
    double L01=0.0;
    double L15=0.0;
    double L04=0.0;
    double L45=0.0;
    double L37=0.0;
    double L23=0.0;
    double L26=0.0;
    double L67=0.0;
    
    double b0,b1,b2,b3;
    double v0,v1,v2,v3,vhex;
    double H12,H47,H30,H56;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i,j+1)-1;
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        L01 = sqrt((P[0*3+0]-P[1*3+0])*(P[0*3+0]-P[1*3+0])+
                   (P[0*3+1]-P[1*3+1])*(P[0*3+1]-P[1*3+1])+
                   (P[0*3+2]-P[1*3+2])*(P[0*3+2]-P[1*3+2]));
        
        L15 = sqrt((P[1*3+0]-P[5*3+0])*(P[1*3+0]-P[5*3+0])+
                   (P[1*3+1]-P[5*3+1])*(P[1*3+1]-P[5*3+1])+
                   (P[1*3+2]-P[5*3+2])*(P[1*3+2]-P[5*3+2]));
        
        H12 = sqrt((P[1*3+0]-P[2*3+0])*(P[1*3+0]-P[2*3+0])+
                   (P[1*3+1]-P[2*3+1])*(P[1*3+1]-P[2*3+1])+
                   (P[1*3+2]-P[2*3+2])*(P[1*3+2]-P[2*3+2]));
        
        
        b0 = 0.5*L01*L15;
        v0 = 1.0/3.0*b0*H12;
        
        L04 = sqrt((P[0*3+0]-P[4*3+0])*(P[0*3+0]-P[4*3+0])+
                   (P[0*3+1]-P[4*3+1])*(P[0*3+1]-P[4*3+1])+
                   (P[0*3+2]-P[4*3+2])*(P[0*3+2]-P[4*3+2]));
        
        L45 = sqrt((P[4*3+0]-P[5*3+0])*(P[4*3+0]-P[5*3+0])+
                   (P[4*3+1]-P[5*3+1])*(P[4*3+1]-P[5*3+1])+
                   (P[4*3+2]-P[5*3+2])*(P[4*3+2]-P[5*3+2]));
        
        H47 = sqrt((P[4*3+0]-P[7*3+0])*(P[4*3+0]-P[7*3+0])+
                   (P[4*3+1]-P[7*3+1])*(P[4*3+1]-P[7*3+1])+
                   (P[4*3+2]-P[7*3+2])*(P[4*3+2]-P[7*3+2]));
        
        b1 = 0.5*L04*L45;
        v1 = 1.0/3.0*b1*H47;
        
        L37 = sqrt((P[3*3+0]-P[7*3+0])*(P[3*3+0]-P[7*3+0])+
                   (P[3*3+1]-P[7*3+1])*(P[3*3+1]-P[7*3+1])+
                   (P[3*3+2]-P[7*3+2])*(P[3*3+2]-P[7*3+2]));
        
        L23 = sqrt((P[2*3+0]-P[3*3+0])*(P[2*3+0]-P[3*3+0])+
                   (P[2*3+1]-P[3*3+1])*(P[2*3+1]-P[3*3+1])+
                   (P[2*3+2]-P[3*3+2])*(P[2*3+2]-P[3*3+2]));
        
        H30 = sqrt((P[3*3+0]-P[0*3+0])*(P[3*3+0]-P[0*3+0])+
                   (P[3*3+1]-P[0*3+1])*(P[3*3+1]-P[0*3+1])+
                   (P[3*3+2]-P[0*3+2])*(P[3*3+2]-P[0*3+2]));
        
        b2 = 0.5*L37*L23;
        v2 = 1.0/3.0*b2*H30;
        
        L26 = sqrt((P[2*3+0]-P[6*3+0])*(P[2*3+0]-P[6*3+0])+
                   (P[2*3+1]-P[6*3+1])*(P[2*3+1]-P[6*3+1])+
                   (P[2*3+2]-P[6*3+2])*(P[2*3+2]-P[6*3+2]));
        
        L67 = sqrt((P[6*3+0]-P[7*3+0])*(P[6*3+0]-P[7*3+0])+
                   (P[6*3+1]-P[7*3+1])*(P[6*3+1]-P[7*3+1])+
                   (P[6*3+2]-P[7*3+2])*(P[6*3+2]-P[7*3+2]));
        
        H56 = sqrt((P[5*3+0]-P[6*3+0])*(P[5*3+0]-P[6*3+0])+
                   (P[5*3+1]-P[6*3+1])*(P[5*3+1]-P[6*3+1])+
                   (P[5*3+2]-P[6*3+2])*(P[5*3+2]-P[6*3+2]));
        
        b3 = 0.5*L26*L67;
        v3 = 1.0/3.0*b3*H56;
        vhex = v0+v1+v2+v3;
    
        vol_cells[i] = vhex;
        
    }
    return vol_cells;
}
*/
