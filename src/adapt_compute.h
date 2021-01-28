#include "adapt_datastruct.h"
#include "adapt_partition.h"
#include "adapt_datatype.h"
#include "adapt_math.h"
#include "adapt.h"

#ifndef ADAPT_COMPUTE_H
#define ADAPT_COMPUTE_H


void NegateVec3D(Vec3D* a);

double DotVec3D(Vec3D* a, Vec3D* b);

Array<double>* MatInv(Array<double>* A);

Array<double>* MatMul(Array<double>* A, Array<double>* B);

double ComputeDetJac(double *P0,double *P1,double *P2,double *P3);

/* 
This needs to be done when linking C programs to Fortran. The reason is name mangling
AKA name decoration
*/
extern "C" {
double ComputeJ(double*P, int ElType);
}

inline double ComputeEdgeLength(Vert* v0, Vert* v1);

double ComputeVolumeHexCell(double *P);

double ComputeDeterminantJ_tet_v2(double*P);
double* ComputeJAtCenter_tet_v2(double*P);

// This function outputs J as an array of 9 values where the matrix is defined as:

/*
 Jac = [J[0], J[1], J[2]
        J[3], J[4], J[5]
        J[6], J[7], J[8]]
*/
// J is computed using the 8-point isoparametric mapping for a hex. The 8-point rule should be sufficient since everything is linear anyways.
Vert* ComputeCenterCoord(double*P, int np);

extern "C" {
double* ComputeJAtCenter(double*P, int np);
}

double ComputeDeterminantJ_tet(double*P);

double ComputeDeterminantJ(double*P, int np);

double ComputeDeterminantJ(double*P, int np);

Array<double>* ComputeDeterminantofJacobian(ParArray<int>* ien, Array<double>* xcn);

double* ComputeVolumeCells(Array<double>* xcn, Array<int>* ien, MPI_Comm comm);

double ComputeTetVolume(double *P);

double* ComputeVolumeCellsReducedToVerts(Array<double>* xcn, Array<int>* ien);

void UnitTestJacobian();

Array<double>* ComputeMetric(Partition* Pa,
                             std::vector<double> d2udx2_v,
                             std::vector<double> d2udxy_v,
                             std::vector<double> d2udxz_v,
                             std::vector<double> d2udyx_v,
                             std::vector<double> d2udy2_v,
                             std::vector<double> d2udyz_v,
                             std::vector<double> d2udzx_v,
                             std::vector<double> d2udzy_v,
                             std::vector<double> d2udz2_v, MPI_Comm comm);

Array<double>* ComputeFaceValues(Partition* P, Array<double>* U, MPI_Comm comm);

Array<double>* ComputeVolumes(Partition* Pa);

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
#endif
