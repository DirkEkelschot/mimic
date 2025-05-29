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

        FaceSetPointer getAllSharedFaceMap();
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

        std::map<int,std::vector<int> >         m_Face2Elem;

        FaceSetPointer                          m_AllSharedFaceSetPointer;

        int vloc; // running total number of vertices.
    };


    #endif
