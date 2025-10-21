#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate_lite.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"
#include "adapt_distri_parstate.h"
#include "adapt_prismtetratrace.h"
#include "adapt_io.h"
#include "adapt_elements.h"

#ifndef ADAPT_PARTOBJECT_LITE_H
#define ADAPT_PARTOBJECT_LITE_H

class PartObjectLite{
        public:
        PartObjectLite(){};
        PartObjectLite(std::map<int,std::vector<int> > Elem2Vert_uniform,
            std::map<int,std::vector<double> > vertices, 
            std::map<int,int> eltype_map,
            std::vector<int> elTypes,
            FaceSetPointer allbcFaces,
            int nE,
            int nV,
            MPI_Comm comm);

        ~PartObjectLite();


        void DeterminePartitionLayout( std::map<int,std::vector<int> > Elem2Vert_uniform,
            std::vector<int> elTypes,
            MPI_Comm comm);


        void DetermineElement2ProcMap(std::map<int,std::vector<int> >   Elem2Vert_uniform, 
                std::map<int,std::vector<double> > vertices_i,
                std::map<int,int>                   Elem2Type_uniform, 
                MPI_Comm comm);

        std::map<int,std::vector<int> > CommunicateAdjacencyInfoLocalPartition(MPI_Comm comm);

        std::set<int> getLocalElemSet();
       
        std::map<int,std::vector<int> > getElem2VertMap();
        
        std::map<int,std::vector<double> > getLocalVertsMap();

        void AddStateVecForAdjacentElements(std::map<int,std::vector<double> > &U, int nvar, MPI_Comm comm);

        std::map<int, std::vector<int> > getAdjacentElementLayer(std::map<int,std::vector<double> > xcn, MPI_Comm comm);

        void GenerateFace2ElementMap(std::map<int,std::vector<int> > Face2Elem_i,  MPI_Comm comm);

        void ComputeFaceMap(MPI_Comm comm, FaceSetPointer allbcFaces);

        void buildExtendedAdjacencyData(MPI_Comm comm, 
                                        std::map<int,std::vector<double> > ghosts, int max_layers);

        FaceSetPointer getAllSharedAndInterFaceFaceMap();

        FaceSetPointer getOwnedSharedAndInterFaceFaceMap();
        std::map<int,std::vector<double> > GatherSharedVertCoordsOnRoot(MPI_Comm comm);
        FaceSetPointer getActualSharedFaceMap();
        std::map<int,int> getActualSharedVerts_Global2LocalMap();
        std::map<int,int> getActualSharedVerts_Local2GlobalMap();
        std::set<int> getOwnedSharedVertsMap();
        std::set<int> getOwnedNonSharedVertsMap();
        FaceSetPointer getExternalFacesForRankFaceMap();
        
        FaceSetPointer getOwnedInteriorFaceFaceMap();
        FaceSetPointer getOwnedBoundaryFaceFaceMap();

        FaceSetPointer getExternalInterFaceFaceMap();
        FaceSetPointer getSharedFacesForRankMap();
        
        

        std::vector<int> getownedSharedFacesOffsets();
        std::vector<int> getownedSharedFacesNlocs();

        std::vector<int> getownedInteriorFacesOffsets();
        std::vector<int> getownedInteriorFacesNlocs();

        std::map<int,int> getLocalVert2GlobalVert();
        std::map<int,int> getGlobalVert2LocalVert();
            
        std::map<int,std::vector<int> > getElem2LocalVertMap();

        void BuildPMMGCommunicationData(MPI_Comm comm, FaceSetPointer allbcFaces);

        int** getParMMGCommFace2GlobalVertMap();
        int** getParMMGCommFace2LocalVertMap();
        int* getParMMGCommColorFace();
        int* getParMMGCommNFacesPerColor();
        int getParMMGNComm();
        std::vector<int> getFace4ParMMG();
        std::map<int,int> getLocalSharedFace2GlobalSharedFace();
        std::map<int,int> getGlobalSharedFace2LocalSharedFace();
        std::map<int,std::vector<int> > getSharedFaceMap();
        std::map<int,std::vector<double> > RedistributeVertexDataForPartition(MPI_Comm comm, int m_Nv_glob, std::map<int, std::vector<double> > t_hessian);
        
        private:
        
        int m_rank;
        int m_size;
        int m_nElemGlobalPart;
        int m_Ne_glob;
        int m_Nf_glob;
        int m_Nv_glob;

        std::vector<int>                        m_part;
        std::map<int,int>                       m_partGlobalRoot;
        std::map<int,int>                       m_partMap;

        std::set<int>                           m_vertIDs_on_rank;
        std::map<int,int>                       m_elem2type_on_rank;

        std::map<int,std::vector<double> >      m_LocalVertsMap;
        std::map<int,int>                       m_LocalV2GlobalV;
        std::map<int,int>                       m_GlobalV2LocalV;
        
        std::set<int>                           m_ElemSet;

        std::map<int,std::vector<double> >      Elem2Data_uniform;
        std::map<int,std::vector<int> >         m_GlobalVert2Elem;

        std::map<int,std::vector<int> >         m_Elem2Face;
        std::map<int,std::vector<int> >         m_Elem2Elem;
        std::map<int,std::vector<int> >         m_Elem2Data;
        std::map<int,std::vector<int> >         m_Elem2Vert;
        std::map<int,std::vector<double> >      m_Elem2Centroid;
        std::map<int,int>                       m_Elem2Rank;
        std::map<int,std::vector<int> >         m_Elem2LocalVert;
        std::map<int,std::vector<int> >         m_Face2Elem;
        FaceSetPointer                          m_SharedFacesForRank;
        FaceSetPointer                          m_AllSharedFaceSetPointer;
        FaceSetPointer                          m_OwnedSharedFaceSetPointer;
        FaceSetPointer                          m_NotOwnedSharedFaceSetPointer;
        FaceSetPointer                          m_ActualSharedFaceSetPointer;
        FaceSetPointer                          m_OwnedInteriorFaceSetPointer;
        FaceSetPointer                          m_ExternalFaceForRankFaceSetPointer;
        FaceSetPointer                          m_OwnedBoundaryFaceSetPointer;
        FaceSetPointer                          m_ExternalInterFaceFace;
        FaceSetPointer                          m_OverallFaceMapOnRank;
        std::map<int,int>                       m_ActualSharedVerts_g2l;
        std::map<int,int>                       m_ActualSharedVerts_l2g;
        std::set<int>                           m_OwnedSharedVerts;
        std::set<int>                           m_OwnedNonSharedVerts;
        std::vector<int>                        m_ownedSharedFacesOffsets;
        std::vector<int>                        m_NownedSharedFaces;
        std::vector<int>                        m_ownedInteriorFacesOffsets;
        std::vector<int>                        m_NownedInteriorFaces;
        std::vector<int>                        m_faces4parmmg;
        std::map<int,std::vector<int> >         m_ColorsFaces;
        std::map<int,std::vector<int> >         m_Elem2AdjElem;
        int vloc; // running total number of vertices.

        int m_ncomm;
        int** m_ifc_tria_glo;
        int** m_ifc_tria_loc;
        int* m_color_face;
        int* m_ntifc;

        std::map<int,int> m_locShF2globShF;
        std::map<int,int> m_globShF2locShF;

        std::map<int,std::vector<int> > m_sharedFaceVerts;
    };


    #endif
