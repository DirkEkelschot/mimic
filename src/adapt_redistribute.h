#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"
#include "adapt_distri_parstate.h"


#ifndef ADAPT_REDISTRIBUTE_H
#define ADAPT_REDISTRIBUTE_H

class RedistributePartitionObject{
    public:
        RedistributePartitionObject(){};
    
        RedistributePartitionObject(US3D* us3d,
                          std::map<int,std::vector<int> > tetras,
						  std::map<int,std::vector<int> > iferank_part_map,
                          std::map<int,std::vector<int> > ief_part_map,
                          std::map<int,std::vector<int> > ifn_part_map,
                          std::map<int,std::vector<int> > ife_part_map,
                          std::map<int,int> if_ref_part_map,
                          std::map<int,std::vector<int> > ushell,
                          std::map<int,Array<double>* > M_vmap,
                          MPI_Comm comm);
    
        ~RedistributePartitionObject();
    
        void RebasePartitionObject(std::map<int,std::vector<int> > tetras,
								   std::map<int,std::vector<int> > iferank_part_map,
                                   std::map<int,std::vector<int> > ief_part_map,
                                   std::map<int,std::vector<int> > ifn_part_map,
                                   std::map<int,std::vector<int> > ife_part_map,
                                   std::map<int,int> if_ref_part_map,
                                   std::map<int,std::vector<int> > ushell,
                                   std::map<int,Array<double>* > M_vmap);
    
        void GetNewGlobalPartitioningTetrahedraMesh();
    
        int* CommunicatePartitionLayout();
        
        void UpdateTetrahedraOnPartition(int nglob,
                                         int nElexpctd,
                                         ParArray<double>* xcn,
                                         ParallelState* xcn_pstate,
                                         ParallelState* ife_pstate);
    
        void GetFace2NodeRedistributedMesh(ParArray<int>* ife, ParArray<int>* if_ref, int ncol, ParallelState* ife_pstate, int nGlob, std::map<int,std::vector<int> > ushell, MPI_Comm comm);
    
        void GetFace2RankTetrahedraMesh();
    
    
        std::map<int,Array<double>* > GetVert2MetricMap();
        Array<int>* GetElement2NodeMap();
        std::map<int,std::vector<int> > GetFace2NodeMap();
        std::map<int,std::vector<int> > GetFace2ElementMap();
        std::vector<Vert*> GetLocalVertices();
        int** GetFace2LocalNode();
        int** GetFace2GlobalNode();
        int GetNBoundaryFaces();
        std::vector<int> GetFaces4ParMMG();
        std::map<int,int> GetGlobalVert2LocalVertMap();
        std::map<int,int> GetLocalVert2GlobalVertMap();
        int* GetColorFace();
        int* GetNFacesPerColor();
        int GetNcomm();
        std::map<int,int> GetLocalSharedFace2GlobalSharedFace();
        std::map<int,int> GetFace2RefMap();
        std::map<int,int> GetShellTet2HybFaceMap();
        std::map<int,std::set<int> > GetShellFace2VertRefMap();
        std::map<int,int > GetShellVert2RefMap();
        std::map<int,int> GetTetF2HybFMap();
        std::map<int,int> GetTag2TetVertMap();
        std::map<int,int> GetTet2TagVertMap();
        std::map<int,int> GetShellVert2RefMap_Global();
        std::map<int,int> GetVertTag2RefLocalMap();
        std::map<int,int > GetShellVert2RefMap_Local();
        std::map<int,int > GetShellVertTag2RefMap_Global();
        std::map<int,Vert*> GetShellVertCoords2RefMap_Global();
        std::map<int,std::vector<int> > GetBndRef2FaceMap();
    private:
        MPI_Comm mpi_comm;
        int world_rank;
        int world_size;
        Array<int>* ElGids;
        Array<int>* ien_part_tetra;
        Array<int>* ien_part_hybrid;
        Array<int>* ief_part_tetra;
        Array<int>* ief_part_hybrid;
        std::map<int,std::vector<int> > face2node;
        std::map<int,std::vector<int> > face2element;

    
        Array<int>* iefref_part_tetra;
    
        std::vector<Vert*> LocalVerts;
        std::map<int,int> m_locV2globV;
        std::map<int,int> m_globV2locV;
    
        std::vector<int*> LocalFaces;
        std::map<int,int> locF2globF;
        std::map<int,int> globF2locF;
    
        std::map<int,int> hybF2tetF;
        std::map<int,int> tetF2hybF;
    
        std::map<int,int> hybV2tetV;
        std::map<int,int> tetV2hybV;
    
        std::map<int,int> hybE2tetE;
        std::map<int,int> tetE2hybE;
    
        std::map<int,int> m_face2ref;
        std::map<int,std::vector<int> > ref2face;
        std::map<int,Array<double>* > m_M_vmap;

        std::map<int,std::vector<int> > m_Boundary_Ref2Face;
        Array<int>* m_part;
    
    
        std::map<int,std::vector<int> > m_TetraToSendToRanks;
        std::map<int,std::vector<int> > m_HybridToSendToRanks;
        std::map<int,std::vector<int> > m_TetraVertIDToSendToRanks;
        std::map<int,std::vector<int> > m_HybridVertIDToSendToRanks;
        std::map<int,std::vector<int> > m_TetraFaceIDToSendToRanks;
        std::map<int,std::vector<int> > m_TetraFaceRefToSendToRanks;
        std::map<int,std::vector<int> > m_HybridFaceIDToSendToRanks;
        std::map<int,std::vector<int> > m_TetraRank2ReqVerts;
        std::map<int,std::vector<int> > m_HybridRank2ReqVerts;
        std::map<int,std::vector<double> > metricsToSend;
        
        std::vector<int> m_TetraVertsOnRank;
        std::vector<int> m_HybridVertsOnRank;

        int m_nPartBoundFaces;
        std::vector<int> m_faces4parmmg;
        std::map<int,int> m_globShF2locShF;
        std::map<int,int> m_locShF2globShF;
        int** m_ifc_tria_glo;
        int** m_ifc_tria_loc;
        int* m_color_face;
        int* m_ntifc;
        int m_ncomm;
        std::map<int,std::vector<int> > m_ColorsFaces;
        std::map<int,int> m_TetEl2HybEl;
        int m_nnpe;
        std::map<int,int> m_shell_tet2hyb_loc;
        std::map<int,int> m_shellVert2RefMapGlobal;
        std::map<int,int> m_shellVert2RefMapLocal;
        std::map<int,int> m_shellVertTag2RefMapLocal;
        std::map<int,int> m_shellVertTag2RefMapGlobal;
        std::map<int,int> m_shellTagID2TetVID;
        std::map<int,int> m_VertTag2RefLocal;
        std::map<int,int> m_ShellVidConsidered;
        std::map<int,std::vector<int> > m_ShellFace2Vert;
        std::map<int,std::set<int> > m_ShellFace2VertRef;
        std::map<int,Vert*> m_shellVertCoords2Ref;
};

#endif

//TetrahedraMesh* ExtractTetrahedralMesh(Array<int>* part_global,
//                                       std::map<int,std::vector<int> > tetras,
//                                       std::map<int,std::vector<int> > ief_part_map,
//                                       std::map<int,std::vector<int> > ifn_part_map,
//                                       std::map<int,std::vector<int> > ife_part_map,
//                                       std::map<int,std::vector<int> > if_ref_part_map,
//                                       std::set<int> ushell,
//                                       std::map<int,Array<double>* > M_vmap,
//                                       MPI_Comm comm)
