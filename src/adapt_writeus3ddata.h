#include "adapt.h"
#include "adapt_repartition.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"
#include "adapt_operations.h"
#include "adapt_output.h"

void WritePrismsUS3DFormat(MPI_Comm comm, 
                            RepartitionObject* prism_repart,
                            PrismTetraTrace* pttrace,
                            std::map<int,int> tracetagV2globalV,
                            std::vector<int> &ifn_P,
                            std::map<int,std::vector<int> > ranges_id);

void WriteBoundaryDataUS3DFormat(MPI_Comm comm, 
                                RepartitionObject* prism_repart,
                                PrismTetraTrace* pttrace,
                                std::vector<int> ifn_T,
                                std::vector<int> ifn_P,
                                std::map<int,std::vector<int> > bcref2bcface_T,
                                std::map<int,int> unique_trace_verts2refmap,
                                std::map<int,int> tracetagV2globalV,
                                std::map<int,std::vector<int> > BoundaryFaces_T,
                                std::map<int,int> glob2locVid_T,
                                std::vector<std::vector<int> > new_tetrahedra_T,
                                std::vector<std::vector<double> > new_vertices_T,
                                std::map<int,int> LocationSharedVert_T,
                                std::map<int,int> oldglob2newglob_T,
                                std::map<int,int> lh_T_bc,
                                std::map<int,int> new_locE2globE_T,
                                std::vector<std::vector<int> > &bcArrays,
                                std::map<int,int> &bcsizing,
                                std::vector<std::vector<int> > zdefs,
                                std::map<int,int> zone2bcref,
                                std::map<int,char*> zone2name,
                                std::vector<std::vector<char> > znames);