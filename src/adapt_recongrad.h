#include "adapt_partition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_topology.h"

#ifndef ADAPT_RECONGRAD_H
#define ADAPT_RECONGRAD_H


std::map<int,Array<double>* > ComputedUdx_LSQ_Vrt_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, 
                                                       std::map<int,double> Uv, 
                                                       Mesh_Topology* meshTopo, 
                                                       Array<double>* ghost, 
                                                       MPI_Comm comm, int world_rank, int world_size);

std::map<int,Array<double>* >  ComputedUdx_LSQ_US3D(Partition* Pa, std::map<int,Array<double>* > U, Array<double>* ghost, MPI_Comm comm, int rank, int size);

std::map<int,Array<double>* > ComputedUdx_MGG(Partition* Pa, std::map<int,double> U,
                               Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm, int rank, int size);
#endif
