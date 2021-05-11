#include "adapt_partition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_topology.h"

#ifndef ADAPT_RECONGRAD_H
#define ADAPT_RECONGRAD_H

std::map<int,Array<double>* > Py_ComputedUdx_LSQ_US3D(std::vector<Vert* > LocalVs,
                                                      std::map<int,std::vector<int> > gE2lV,
                                                      std::map<int,int> gV2lV,
                                                      std::vector<int> Loc_Elem,
                                                      std::map<int,std::vector<int> > ifn_part_map,
                                                      std::map<int,std::vector<int> > ief_part_map,
                                                      std::map<int,std::vector<int> > iee_part_map,
                                                      std::map<int,std::vector<int> > if_Nv_part_map,
                                                      std::map<int,int> LocElem2Nf,
                                                      std::map<int,int> LocElem2Nv,
                                                      int Nel_glob,
                                                      std::map<int,Array<double>* > UState,
                                                      Array<double>* ghost, MPI_Comm comm);

std::map<int,Array<double>* > ComputedUdx_LSQ_Vrt_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, std::map<int,double> Uv, Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);

std::map<int,Array<double>* >  ComputedUdx_LSQ_US3D(Partition* Pa, std::map<int,Array<double>* > U, Array<double>* ghost, MPI_Comm comm);

std::map<int,Array<double>* > ComputedUdx_MGG(Partition* Pa, std::map<int,double> U,
                               Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);
#endif
