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

double ComputeEdgeLength(Vert* v0, Vert* v1);

double ComputeVolumeHexCell(double *P);
double ComputeVolumeTetCell(double *P);
double ComputeVolumePrismCell(double *P);

double ComputeDeterminantJ_tet_v2(double*P);
Array<double>* ComputeJAtCenter_tet_v2(double*P);

// This function outputs J as an array of 9 values where the matrix is defined as:

/*
 Jac = [J[0], J[1], J[2]
        J[3], J[4], J[5]
        J[6], J[7], J[8]]
*/
// J is computed using the 8-point isoparametric mapping for a hex. The 8-point rule should be sufficient since everything is linear anyways.
Vert* ComputeCenterCoord(double*P, int np);
Vert* ComputeCentroidCoord(double*P, int np);

extern "C" {
double* ComputeJAtCenter(double*P, int np);
}

double ComputeDeterminantJ_tet(double*P);

double ComputeDeterminantJ(double*P, int np);

double ComputeDeterminantJ(double*P, int np);

Array<double>* ComputeDeterminantofJacobian(ParArray<int>* ien, Array<double>* xcn);

double* ComputeVolumeCells(Array<double>* xcn, Array<int>* ien, MPI_Comm comm, int rank, int size);

double ComputeTetVolume(double *P);

double* ComputeVolumeCellsReducedToVerts(Array<double>* xcn, Array<int>* ien);

void UnitTestJacobian();

void ComputeMetric(Partition* Pa, std::vector<double> metric_inputs, MPI_Comm comm, int rank, int size,
                   std::map<int,Array<double>* > scale_vm,
                   std::map<int,Array<double>* > &Hess_vm,
                   double sumvol, double po);

Array<double>* ComputeFaceValues(Partition* P, Array<double>* U, MPI_Comm comm, int rank, int size);

Array<double>* ComputeVolumes(Partition* Pa);

#endif
