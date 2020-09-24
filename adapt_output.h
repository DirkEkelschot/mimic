#include "adapt_datastruct.h"
//#include "adapt_part_func.h"
#include "adapt_partition.h"
#include "adapt_io.h"
#include "adapt.h"
#include "adapt_array.h"
#include "adapt_topology.h"
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"
#ifndef ADAPT_OUTPUT_H
#define ADAPT_OUTPUT_H

using namespace std;

void OutputMesh_MMG(MMG5_pMesh mmgMesh);

void OutputBoundaryID_MMG(MMG5_pMesh mmgMesh, std::map<int,std::vector<int> > ref2bface, int bndID);

void OutputBoundaryID(Partition* P, Mesh_Topology* meshTopo, US3D* us3d, int bndID);

void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs,MPI_Comm comm);

void OutputPartition(Partition* part, ParArray<int>* ien, Array<double>* H,  MPI_Comm comm);

void OutputCompletePartition(Partition* part, ParArray<int>* ien, Array<double>* H, MPI_Comm comm);

void OutputZone(Partition* part, Array<double>* H, MPI_Comm comm);

void OutputGradient(Partition* parttn, Array<double>* H, ParallelState* pstate, MPI_Comm comm);

void OutputQuantityPartition(Partition_old* pa, Array<double>* Quan, MPI_Comm comm);

//void OutputPartionVolumes(ParArray<int>* ien, Array<double>* xcn_on_root, MPI_Comm comm);

//void OutputPartitionFaces();

void WriteBoundaryDataInSerial3(Array<double>* xcn);

void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts);

void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts);


#endif
