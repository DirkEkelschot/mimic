#include "adapt_partition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_topology.h"

#ifndef ADAPT_RECONGRAD_H
#define ADAPT_RECONGRAD_H

Array<double>* ComputedUdx_LSQ(Partition* P, std::vector<double> U, int Nel, MPI_Comm comm);

Array<double>* ComputedUdx_LSQ_US3D_v1(Partition* P, ParallelState* pstate, ParArray<int>* iee, std::map<int,std::vector<int> > iee_vec,std::map<int,std::vector<int> > ief_vec, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, Array<double>* bound, MPI_Comm comm);

Array<double>* ComputedUdx_LSQ_US3D_v2(Partition* P, ParallelState* pstate, ParArray<int>* iee, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, Array<double>* bound, MPI_Comm comm);

std::map<int,Array<double>* >  ComputedUdx_LSQ_US3D_v3(Partition* Pa, std::map<int,double> U,Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);

Array<double>* ComputedUdx_MGG(Partition* Pa, std::map<int,double> U,
                               Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);
#endif
