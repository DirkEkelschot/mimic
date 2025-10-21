
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

#ifndef ADAPT_PARTOBJECT_H
#define ADAPT_PARTOBJECT_H

class PartObject{
        public:
                PartObject(){};

                PartObject(mesh* meshInput,
                    std::map<int,std::vector<int> > Elem2Vert_i,
                    std::map<int,std::vector<int> > Elem2Face_i,
                    std::map<int,std::vector<int> > Elem2Elem_i,
                    std::map<int,std::vector<double> > Elem2Data_i,
                    std::map<int,int> Elem2Type_i,
                    int nAdjLayer,
                    bool reconstruct_ifn,
                    MPI_Comm comm);

                ~PartObject();

                void GetOptimalDistributionSchedule(std::map<int,std::vector<int> > Elem2Vert_i,
                                                        std::map<int,std::vector<int> > Elem2Face_i,
                                                        std::map<int,std::vector<int> > Elem2Elem_i,
                                                        std::map<int,std::vector<double> > Elem2Data_i,
                                                        std::map<int,int> Elem2Type_i,
                                                        std::map<int,std::vector<int> >& Elem2Vert_uniform,
                                                        std::map<int,std::vector<int> >& Elem2Face_uniform,
                                                        std::map<int,std::vector<int> >& Elem2Elem_uniform,
                                                        std::map<int,std::vector<double> >& Elem2Data_uniform,
                                                        std::map<int,int>& Elem2Type_uniform,
                                                        MPI_Comm comm);
                
        
                void DeterminePartitionLayout( std::map<int,std::vector<int> > Elem2Vert_uniform,
                                                std::vector<int> element2rank,
                                                std::vector<int> elTypes,
                                                MPI_Comm comm);

                void DetermineElement2ProcMap(std::map<int,std::vector<int> >     Elem2Vert_uniform, 
                                                std::map<int,std::vector<int> >     Elem2Face_uniform,
                                                std::map<int,std::vector<int> >     Elem2Elem_uniform,
                                                std::map<int,std::vector<double> >  Elem2Data_uniform,
                                                std::map<int,int>                   Elem2Type_uniform, 
                                                std::map<int,std::vector<double> >  vertices_i,
                                                int Nf_glob,
                                                int Nv_glob, 
                                                MPI_Comm comm);

                void GenerateFace2ElementMap(std::map<int,std::vector<int> > Face2Elem_i,
                                            MPI_Comm comm);

                void GenerateFace2VertexMap(MPI_Comm comm);

                void getFace2VertexPerPartitionMap(std::map<int,std::vector<int> > ifn_read, 
                                               std::map<int,std::vector<int> > if_Nv_read, 
                                               int Nf_glob, 
                                               MPI_Comm comm);

                std::map<int, std::vector<int> > CommunicateAdjacencyInfoLocalPartition(MPI_Comm comm);
                std::map<int, std::vector<int> > CommunicateAdjacencyInfoExtendedPartition(MPI_Comm comm);

                void GenerateTraceMap();

                std::map<int, std::vector<int> > getAdjacentElementLayer(std::map<int,std::vector<double> > xcn, MPI_Comm comm);


                void updateFace2EntityPerPartition(std::map<int,std::vector<int> > adjacent_ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > if_Nv_read,
                                                   MPI_Comm comm);
                
                void buildExtendedAdjacencyDataPerLayer(MPI_Comm comm, 
                                            std::map<int,std::set<int> >& E2A_set,
                                            std::set<int> input_elems,
                                            std::map<int,std::set<int> >& toCheck,
                                            std::set<int>& toCheckset,
                                            std::map<int,std::vector<double> > ghosts);

                void buildExtendedAdjacencyData(MPI_Comm comm, 
                                                std::map<int,std::vector<double> > ghosts);
                        
                void buildExtendedAdjacencyDataReFactor(MPI_Comm comm, 
                                           std::map<int,std::vector<double> > ghosts, int max_layers);

                std::map<int,std::vector<int> > getExtendedAdjacencyDataReFactor();

                std::map<int,int> getPartMap();
                std::map<int,std::vector<int> > getElem2VertMap();
                std::map<int,std::vector<int> > getElem2FaceMap();
                std::map<int,std::vector<int> > getFace2VertMap();
                std::map<int,std::vector<int> > getElem2ElemMap();
                std::map<int,std::vector<int> > getFace2ElemMap();
                std::map<int,int> getLocalVert2GlobalVert();
                std::map<int,std::vector<double> > getElem2DataMap();
                std::map<int,std::vector<int> > getElem2LocalVertMap();
                std::map<int,std::vector<double> > getElem2CentroidData();
                std::map<int,std::vector<int> > getExtendedAdjacencyData();

                std::map<int,int > getElem2RankMap();
                std::set<int> getTraceVertsOnRankMap();
                std::set<int> getTraceFacesOnRankMap();
                std::set<int> getLocalElemSet();
                std::map<int,std::vector<double> > getLocalVertsMap();
                std::map<int,std::vector<double> > getGhostFaceVert();

                std::map<int,std::set<int> > getExtendedAdjacencyData(MPI_Comm comm, std::map<int,std::vector<double> > &GhostVerts);

                void AddStateVecForAdjacentElements(std::map<int,std::vector<double> > &U, int nvar, MPI_Comm comm);

                void SetStateVec(std::map<int,std::vector<double> > U, int nvar);

                std::map<int,std::map<int,double> > GetNode2ElementMap();

                std::map<int,std::vector<double> > ReduceStateVecToVertices(std::map<int,std::map<int,double> > Vert2ElemMap,
                                                std::map<int,std::vector<double> > Umap,
                                                int nvar);

                std::map<int,std::map<int, double> > getElem2ConnectedVertMap();

                std::set<int> GetLocalTraceFacesSet();
                std::set<int> GetLocalTraceVertsSet();
                std::map<int,std::vector<int> > GetLocalTraceFace2LeftRight();
                std::map<int,int> getGlobalElement2Rank();
                std::map<int,int> GetElement2TypeOnRankMap();
        

        private:

                int nElemGlobalPart;
                int Ne_glob;
                int Nf_glob;
                int Nv_glob;
                std::vector<int>                        m_part;
                std::map<int,int>                       m_partGlobalRoot;
                std::map<int,int>                       m_partMap;

                // Fundamental data structures (Element2Entity)
                std::map<int,int> m_Elem2Type;
                std::set<int>                           m_ElemSet;
                std::map<int,std::vector<int> >         m_Elem2Vert;
                std::map<int,std::vector<int> >         m_Elem2Face;
                std::map<int,std::vector<int> >         m_Elem2Elem;
                std::map<int,std::vector<double> >      m_Elem2Data;
                std::map<int,std::vector<int> >         m_Face2Elem;
                std::map<int,std::vector<int> >         m_Face2Vert;
                std::map<int,std::vector<double> >      m_LocalVertsMap;
                std::map<int,int>                       m_Elem2Rank;
                std::map<int,std::vector<int> >         m_Rank2Elem;

                // Trace data structures.
                std::map<int,std::vector<int> >         m_TraceFace2Vert;
                std::map<int,std::vector<int> >         m_TraceFace2Elem;
                std::set<int>                           m_TraceVertsOnRank;
                std::set<int>                           m_TraceFacesOnRank;

                std::set<int>                           m_vertIDs_on_rank;
                std::set<int>                           m_faceIDs_on_rank;
                std::map<int,int>                       m_elem2type_on_rank;
                std::map<int,std::set<int> >            m_Rank2ReqElem;

                std::map<int, std::vector<double> >     m_Elem2Centroid;
                std::map<int, std::vector<int> >        m_globElem2globFaces;
                std::map<int, std::vector<int> >        m_globFace2GlobalElements;

                std::map<int,int>                       m_LocalV2GlobalV;
                std::map<int,int>                       m_GlovalV2LocalV;
                std::map<int,std::vector<int> >         m_GlobalVert2Elem;
                std::map<int,std::vector<int> >         m_Elem2LocalVert;
                std::map<int,std::vector<double> >      m_GhostFaceVert;
                std::map<int,std::vector<int> >         m_Elem2AdjElem;
                std::map<int,std::vector<int> >         m_Elem2AdjElem_ReFactor;
                int vloc; // running total number of vertices.
                int floc; // running total number of faces.
};


#endif