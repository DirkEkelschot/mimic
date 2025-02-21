#include "adapt.h"
#include "adapt_repartition.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"
#include "adapt_operations.h"
#include "adapt_prepareadaption.h"

#ifndef ADAPT_TESTING_H
#define ADAPT_TESTING_H

void RunParMMGandWriteTetraUS3Dformat(MPI_Comm comm, 
                                        PMMG_pParMesh &parmesh,
                                        std::map<int,int> part_global_bl,
                                        int bndIDmax,
                                        Inputs* inputs, 
                                        int nElemsGlob_P, int nVertsGlob_P, 
                                        std::map<int,int>& tracetagV2globalV, 
                                        std::vector<int> &ifn_T,
                                        std::map<int,std::vector<int> > &bcmap,
                                        std::map<int,std::vector<int> > &BoundaryFaces_T,
                                        std::map<int,int> &glob2locVid,
                                        std::map<int,int> &LocationSharedVert_update,
                                        std::vector<double> &vertices_output,
                                        std::vector<std::vector<int> > &new_tetrahedra,
                                        std::map<int,int> &oldglob2newglob,
                                        std::map<int,int> &lh_T_bc,
                                        std::map<int,int> &new_globE2locE,
                                        int &nLocIntVrts,
                                        int &nLocShVrts,
                                        std::map<int,int> tagE2gE_P);

PMMG_pParMesh InitializeParMMGmesh(MPI_Comm comm, 
                                    std::set<int> loc_trace_faces,
                                   std::set<int> loc_trace_verts,
                                   std::map<int,std::vector<int> > loc_trace_face2leftright,
                                   PrepareAdaption* prepare,
                                   std::map<int,std::vector<int> > ranges_id,
                                   int bndIDmax,
                                   std::map<int, std::vector<double> > metric_vmap);

void RunParMMGAndTestPartitioning(MPI_Comm comm, 
                                  PMMG_pParMesh parmesh, 
                                  PrepareAdaption* prepare,
                                  std::map<int,std::vector<int> > ranges_id,
                                  Inputs* inputs);


#endif