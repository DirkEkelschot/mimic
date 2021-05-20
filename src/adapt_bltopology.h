#include "adapt.h"
#include "adapt_array.h"
#include "adapt_datatype.h"
#include "adapt_datastruct.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"
#include "adapt_topology.h"
#include "adapt_output.h"
struct BLShellInfo{
  
    std::map<int,int> ShellFace2BFace;
    std::map<int,int> BFace2ShellFace;
    std::map<std::set<int>,int> ShellTri2FaceID;
    std::map<int,std::vector<int> > ShellFaceID2TriID;
    std::map<int,int> FaceID2TopoType;
    std::map<int,std::map<int,int> > ShellFace2ShellVert2OppositeBoundaryVerts;
    Array<int>* ShellRef;
    std::map<int,std::vector<int> > BLlayers;
    std::set<int> elements_set;
};


BLShellInfo* FindOuterShellBoundaryLayerMesh(int wall_id, int nLayer,
                                             Array<double>* xcn_g, 
                                             Array<int>* ien_g, 
                                             Array<int>* ief_g, 
                                             Array<int>* ife_g, 
                                             Array<int>* ifn_g,
                                             ParallelState* xcn_pstate, 
                                             ParallelState* ien_pstate, 
                                             std::map<int,std::vector<int> > bnd_face_map, 
                                             std::map<int,int> vert_ref_map, 
                                             MPI_Comm comm, int world_rank, int world_size);

Mesh_Topology_BL* ExtractBoundaryLayerMeshFromShell(std::vector<std::vector<int> > u_tris, 
                                                    BLShellInfo* BLshell, 
                                                    int wall_id,
                                                    int nLayer, 
                                                    Array<double>* xcn_g, 
                                                    Array<int>* ien_g, 
                                                    Array<int>* ief_g, 
                                                    Array<int>* ife_g,
                                                    Array<int>* ifn_g, 
                                                    ParallelState* xcn_pstate, 
                                                    ParallelState* ien_pstate, 
                                                    std::map<int,std::vector<int> > bnd_face_map, 
                                                    std::map<std::set<int>,int> tria_ref_map, 
                                                    std::map<std::set<int>,int> quad_ref_map,  
                                                    MPI_Comm comm, int world_rank, int world_size);

