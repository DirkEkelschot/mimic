#include "adapt_partition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_topology.h"

#ifndef ADAPT_RECONGRAD_H
#define ADAPT_RECONGRAD_H

std::map<int,Array<double>* > ComputedUdx_LSQ_HO_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, std::map<int,double> gbMap, MPI_Comm comm);

std::map<int,Array<double>* > ComputedUdx_LSQ_Vrt_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, std::map<int,Array<double>* > Uv, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm);

std::map<int,Array<double>* >  ComputedUdx_LSQ_US3D(Partition* Pa, std::map<int,Array<double>* > U, std::map<int,double> gbMap, MPI_Comm comm);

std::map<int,Array<double>* > ComputedUdx_MGG(Partition* Pa, std::map<int,Array<double>*> U,
                               Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm);
#endif
