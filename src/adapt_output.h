#include "adapt_datastruct.h"
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

void OutputBoundaryLayerPrisms(Array<double>* xcn_g, Mesh_Topology_BL* BLmesh, MPI_Comm comm,string fname);

void OutputMesh_MMG(MMG5_pMesh mmgMesh, int offset, int Nel, string fname);

void OutputMesh_MMG_Slice(MMG5_pMesh mmgMesh, int offset, int Nel, string fname);

void OutputBoundaryID_MMG(MMG5_pMesh mmgMesh, std::map<int,std::vector<int> > ref2bface, int bndID);

void OutputBoundaryID(Partition* P, Mesh_Topology* meshTopo, US3D* us3d, int bndID);

void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs);

void OutputPartition(Partition* part, ParArray<int>* ien, Array<double>* H,  MPI_Comm comm, int rank, int size);

void OutputCompletePartition(Partition* part, ParArray<int>* ien, Array<double>* H, MPI_Comm comm, int rank, int size);

void OutputZone(Partition* part, Array<double>* H, MPI_Comm comm, int rank, int size);

void OutputGradient(Partition* parttn, Array<double>* H, ParallelState* pstate, MPI_Comm comm, int rank, int size);

void OutputQuantityPartition(Partition_old* pa, Array<double>* Quan, MPI_Comm comm, int rank, int size);

void OutputBLElements(Partition* part, Mesh_Topology_BL* mesh_topology_bl,  MPI_Comm comm,string fname);

void OutputBLElementsOnRoot(Array<double>* xcn_root, Array<int>* ien_root, std::vector<int> elements,  MPI_Comm comm, int rank, int size, string fname);

void WriteBoundaryDataInSerial3(Array<double>* xcn);

void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts);

void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts);


#endif
