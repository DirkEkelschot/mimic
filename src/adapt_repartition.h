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

#ifndef ADAPT_REPARTITION_H
#define ADAPT_REPARTITION_H

class RepartitionObject{
        public:
                RepartitionObject(){};
                RepartitionObject(mesh* meshInput,
                                  std::map<int,std::vector<int> > elements2verts,
                                  std::map<int,std::vector<int> > elements2faces,
                                  std::map<int,std::vector<int> > elements2elements,
                                  std::map<int,int> element2type,
                                  std::map<int,std::vector<double> > data,
                                  int nAdjLayer,
                                  bool reconstruct_ifn,
                                  MPI_Comm comm);
                ~RepartitionObject();

                void GetSharedTraces(std::map<int,std::vector<int> > elements,
                                        std::map<int,std::vector<int> > ife,
                                        std::map<int,std::vector<int> > if_ref,
                                        std::vector<int> element2rank,
                                        MPI_Comm comm);

                void GetOptimalDistributionSchedule(std::map<int,std::vector<int> > elements2verts,
                                                    std::map<int,std::vector<int> > elements2faces,
                                                    std::map<int,std::vector<int> > elements2elements,
                                                    std::map<int,int> element2type,
                                                    std::map<int,std::vector<double> > data,
                                                    MPI_Comm comm, 
                                                    std::map<int,int>& element2type_update,
                                                    std::map<int,std::vector<int> >& elements2verts_update, 
                                                    std::map<int,std::vector<int> >& elements2faces_update, 
                                                    std::map<int,std::vector<int> >& elements2elements_update,
                                                    std::map<int,std::vector<double> >& elements2data_update);

                void DeterminePartitionLayout(std::map<int,std::vector<int> > elements, std::vector<int> element2rank, std::vector<int> elTypes, MPI_Comm comm);

                
                void DetermineElement2ProcMap(std::map<int,std::vector<int> >     ien, 
                                                 std::map<int,std::vector<int> >   ief,
                                                 std::map<int,std::vector<int> >   iee,  
                                                 std::map<int,std::vector<double> > data,
                                                 std::map<int,int> element2type,
                                                 std::map<int,std::vector<double> >   xcn,
                                                 int Nf_glob,
                                                 int Nv_glob, 
                                                 MPI_Comm comm,
                                                 std::map<int,std::vector<int> >& elements2verts_update, 
                                                 std::map<int,std::vector<int> >& elements2faces_update, 
                                                 std::map<int,std::vector<int> >& elements2elements_update,
                                                 std::map<int,std::vector<double> >& elements2data_update);

                void getFace2EntityPerPartition(std::map<int,std::vector<int> > ief,
                                                std::map<int,std::vector<int> > ife,  
                                                int Nf_glob,
                                                std::map<int,std::vector<int> > &ife_loc,
                                                MPI_Comm comm);

                void getFace2EntityPerPartitionVector(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   std::map<int,std::vector<int> > ifNv_read, 
                                                   std::map<int,std::vector<int> > ifn_read,  
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   std::map<int,std::vector<int> > &ifNv_loc,
                                                   std::map<int,std::vector<int> > &ifn_loc,
                                                   MPI_Comm comm);

                void getFace2EntityPerPartitionRef(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   MPI_Comm comm);

                void getFace2VertexPerPartitionMap(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ifn_read, 
                                                   std::map<int,std::vector<int> > if_Nv_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   MPI_Comm comm);
                
                void getFace2EntityPerPartitionMap(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   std::map<int,std::vector<int> > ifNv_read, 
                                                   std::map<int,std::vector<int> > ifn_read,  
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   std::map<int,std::vector<int> > &ifn_loc,
                                                   MPI_Comm comm);

                void updateFace2EntityPerPartition(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > if_Nv_read,  
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
                                                   MPI_Comm comm);

                
                std::map<int,std::vector<double> > getElement2DataMap();

                std::map<int, std::vector<int> > getGlobalElement2LocalVertMap();

                std::vector<std::vector<double> > getLocalVerts();

                std::map<int, std::vector<double> > getLocalVertsMap();

                std::map<int,int> getGlobalElement2Rank();

                std::map<int,std::vector<int> > getElement2VertexMap();

                std::map<int,std::vector<int> > getElement2ElementMap();

                std::map<int,std::vector<int> > getElement2FacesMap();

                std::map<int, std::vector<int> > getFace2VertexMap();

                std::map<int,std::vector<int> > getFace2ElementMap();

                std::map<int,std::vector<int> > getFace2RefMap();

                std::map<int, std::vector<int> > getFace2NVertexMap();

                std::map<int,int> getTag2ElementTrace();

                std::vector<int> getLocElem();
                
                void buildUpdatedVertexAndFaceNumbering(MPI_Comm comm, 
                                                             std::map<int,std::vector<int> > ranges_id,
                                                             std::map<int,std::vector<int> > ranges_ref);

                std::map<int,int> getVertTag2LocalVert();
                std::map<int,int> getLocalVert2VertTag();

                std::map<int,int> getLocalElement2ElementTag();
                std::map<int,int> getElementTag2LocalElement();
                std::map<int,int> getGlobalElement2ElementTag();
                std::map<int,int> getElementTag2GlobalElement();
                std::map<int,int> getGlob2TagElementID();
                std::map<int,int> getBoundaryFaces();
                std::map<int,int> getBoundaryFaces2Ref();
                std::map<int,int> getUpdatedGlobal2LocalFMap();
                std::map<int,int> getUpdatedLocal2GlobalFMap();
                std::map<int,int> getUpdatedGlobal2LocalVMap();
                std::map<int,int> getUpdatedLocal2GlobalVMap();
                std::map<int,int> getUpdatedTag2GlobalVMap();
                std::map<int,int> getUpdatedGlobal2TagVMap();
                

                void buildCommunicationMaps(MPI_Comm comm);
                int** getParMMGCommFace2GlobalVertMap();
                int** getParMMGCommFace2LocalVertMap();
                int* getParMMGCommColorFace();
                int* getParMMGCommNFacesPerColor();
                int getParMMGNComm();
                
                std::vector<int> getFace4ParMMG();
                std::map<int,int> getGlobal2TagFMap();
                void buildInteriorSharedAndBoundaryFaceMaps(MPI_Comm comm, 
                                                            std::map<int,std::vector<int> > ranges_id,
                                                            std::map<int,std::vector<int> > ranges_ref);
                
                std::map<int, std::vector<int> > getSharedFaceMap();
                std::map<int, std::vector<int> > getInteriorFaceMap();
                std::map<int, std::vector<int> > getBoundaryFaceMap();
                std::map<int,int> getTagF2LocFID();
                // std::map<int,int> getUpdatedLocal2GlobalVertexMap();
                // std::map<int,int> getUpdatedGlobal2LocalVertexMap();
                std::map<int,int> getNonSharedVertsOwned();
                std::map<int,int> getSharedVertsOwned();
                std::map<int,int> getSharedVertsNotOwned();
                std::map<int,int> GetLeftHandFaceElementMap();
                std::map<int,int> GetRightHandFaceElementMap();
                std::map<int,int> GetElement2TypeOnRankMap();
                std::map<int,std::vector<int> > getZone2boundaryFaceID();
                void buildParMMGCommunicationMaps(MPI_Comm comm);

                std::set<int> GetLocalTraceFacesSet();
                std::set<int> GetLocalTraceVertSet();
                std::map<int,std::vector<int> > GetLocalTraceFace2VertMap();
                std::map<int,std::vector<int> > GetLocalTraceFace2LeftRight();

                void ReconstructFace2ElementMap(std::map<int,std::vector<int> > ief,
                                                std::map<int,std::vector<int> > ife_read,
                                                std::map<int,std::vector<int> > &face2elements_update_new,
                                                MPI_Comm comm);

                void ReconstructFace2VertexMap(std::map<int,std::vector<int> > elements2faces_update,
                                  std::map<int,std::vector<int> > &face2verts_update_new,
                                  MPI_Comm comm);


                std::map<int,int> GetLocalSharedFace2GlobalSharedFace();
                // std::map<int,int> GetGlobalSharedFace2LocalSharedFace();

                std::map<int,std::vector<int> > getFaceTag2VertTagMap();
                void buildSharedVertexMap(MPI_Comm comm, PrismTetraTrace* trace);

                std::map<int,std::vector<int> > getAdjacentElementLayer(std::map<int,std::vector<int> > element2verts,
                                             std::map<int,std::vector<int> > element2faces,
                                             std::map<int,std::vector<int> > element2element,
                                             std::map<int,std::vector<double> > xcn, 
                                             std::map<int,std::vector<double> > U, 
                                             int Ne_glob,
                                             int Nf_glob,
                                             int Nv_glob, 
                                             MPI_Comm comm,
                                             std::map<int,std::vector<int> >& elements2verts_update,
                                             std::map<int,std::vector<int> >& elements2faces_update,
                                             std::map<int,std::vector<int> >& elements2elements_update,
                                             std::map<int,std::vector<double> >& elements2data_update);
                
                
                std::map<int,std::set<int> > GetNode2ElementMap();


                void AddStateVecForAdjacentElements(std::map<int,std::vector<double> > &U, int nvar, MPI_Comm comm);
                void SetStateVec(std::map<int,std::vector<double> > U, int nvar);

        private:
                std::vector<int> part;
                std::map<int,int> part_global;
                std::map<int,int> part_map;

                int eloc;
                int vloc;
                int floc;

                int Ne_glob;
                int Nf_glob;
                int Nv_glob;


                std::set<int> m_loc_trace_faces;
                std::set<int> m_loc_trace_verts;
                std::map<int,std::vector<int> > m_loc_trace_face2leftright;
                std::map<int,std::vector<int> > m_loc_trace_face2vertmap;

                std::map<int,std::set<int> > node2elem_map;

                std::map<int, std::vector<int> > globElem2globVerts;
                std::map<int, std::vector<int> > globElem2locVerts;

                std::vector<std::vector<int> > lelement2lvertex;

                std::map<int,int> lface2gface;
                std::map<int,int> gface2lface;

                std::set<int> unique_vertIDs_on_rank_set;
                std::vector<int> unique_verts_on_rank_vec;
                std::set<int> unique_faceIDs_on_rank_set;

                std::vector<int> loc_r_nv_elem;
                std::vector<int> loc_r_nf_elem;
                std::set<int> elem_set;
                std::map<int,int> elem_map;
                std::set<int> loc_r_elem_set;
                std::map<int,std::vector<double> > loc_data;
                std::map<int, std::vector<double> > LocalVertsMap;
                int itel_locadj;
                int nLoc_Elem;
                int nLocAndAdj_Elem;
                int nLoc_Verts;

                // Book keeping maps
                std::map<int,std::vector<int> > globElem2localFaces;
                std::map<int,std::vector<int> > globElem2globFaces;
                std::map<int,std::vector<int> > globFace2GlobalElements;
                std::vector<int> Loc_Elem;
                std::vector<int> Loc_Elem_Nv;
                std::vector<int> Loc_Elem_Nf;
                std::map<int,int> locPartV2globV;
                std::map<int,int> globPartV2locV;
                std::set<int> Loc_Elem_Set;
                std::map<int,std::vector<int> > o_zone2bcface;
                std::map<int,int> elem2type_on_rank;
                std::map<int,int> LocElem2Nv;
                std::map<int,int> LocElem2Nf;
                std::map<int,int> LocElem2Ndata;
                std::vector<int> LocAndAdj_Elem;
                std::vector<int> LocAndAdj_Elem_Nv;
                std::vector<int> LocAndAdj_Elem_Nf;
                std::map<int,int> LocalElement2GlobalElement;
                std::map<int,int> GlobalElement2LocalElement;
                std::map<int, std::vector<int> > globVerts2globElem;

                // Core maps between Elements/Faces/Vertices
                std::map<int, std::vector<int> > elements2verts_update;
                std::map<int, std::vector<int> > elements2faces_update;
                std::map<int, std::vector<int> > elements2elements_update;
                std::map<int, std::vector<double> > elements2data_update;
                std::map<int, std::vector<int> > face2verts_update;
                std::map<int, std::vector<int> > face2Nverts_update;
                std::map<int,std::vector<int> > face2elements_update;
                std::map<int,std::vector<int> > face2reference_update;
                std::map<int,int> o_lvertex2gvertex_part;
                std::map<int,int> o_gvertex2lvertex_part;
                std::map<int,int> o_lvertex2gvertex;
                std::map<int,int> o_gvertex2lvertex;

                //buildInteriorSharedAndBoundaryFacesMaps()
                std::map<int, std::vector<int> > o_sharedFace2Nodes;
                std::map<int, std::vector<int> > o_boundaryFace2Nodes;
                std::map<int, std::vector<int> > o_interiorFace2Nodes;
                std::map<int,int> o_SharedVertsOwned;
                std::map<int,int> o_NonSharedVertsOwned;
                std::map<int,int> o_SharedVertsNotOwned;
                std::map<int,int> o_sharedFace2RankMap;
                std::map<int,int> o_sharedVertexMapUpdatedGlobalID;
                std::map<int,int> o_sharedFaceMapUpdatedGlobalID;
                std::map<int,int> o_sharedVertex2RankMap;   
                std::map<int,int> o_boundaryFaces;
                std::map<int,int> o_boundaryFaces2Ref;
                std::map<int,std::vector<int> > o_SharedFace2Rank;

                std::map<int,int> tag2element_trace;
                std::map<int,std::vector<int> > ref2bcface;
                std::map<int,std::vector<int> > zone2bcface;                

                //buildInteriorSharedAndBoundaryFaceMaps()
                std::map<int,int> o_lhp;
                std::map<int,int> o_rhp;
                std::map<int,int> o_tagE2gE;
                std::map<int,int> o_gE2tagE;
                std::map<int,int> o_tagE2lE;
                std::map<int,int> o_lE2tagE;
                std::map<int,int> o_loc2globF;
                std::map<int,int> o_glob2locF;
                std::map<int,int> o_loc2globV;
                std::map<int,int> o_glob2locV;
                std::map<int,std::vector<int> > o_colorRh;
                std::map<int,std::vector<int> > o_colorFh;
                std::map<int,int> o_tagF2locFID;

                //buildParMMGCommunicationMaps()
                int** o_ifc_tria_glo;
                int** o_ifc_tria_loc;
                int* o_color_face;
                int* o_ntifc;
                int o_ncomm;
                std::map<int,std::vector<int> > o_ColorsFaces;

                //UpdateGlobalIDs()
                std::map<int,std::vector<int> > o_element2verts_global;
                std::map<int,std::vector<int> > o_element2faces_global;
                std::map<int,std::vector<int> > o_face2elements_global;
                std::map<int,std::vector<int> > o_face2verts_global;
                std::map<int,int> o_tag2globV;
                std::map<int,int> o_tag2globF;
                std::map<int,int> o_glob2tagF;
                std::map<int,int> o_glob2tagV;
                std::map<int,int> o_globShF2locShF;

                //GetFace2RankMesh()
                std::vector<int> o_faces4parmmg;

                std::map<int,int> o_locShF2globShF;

                


};      

#endif


