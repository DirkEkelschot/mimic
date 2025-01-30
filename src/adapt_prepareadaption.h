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
#include "adapt_partobject.h"

#ifndef ADAPT_PREPAREADAPTION_H
#define ADAPT_PREPAREADAPTION_H

class PrepareAdaption{
        public:

            PrepareAdaption(){};

            ~PrepareAdaption();

            PrepareAdaption(PartObject* part, 
                            MPI_Comm comm,
                            std::map<int,std::vector<double> >&&   LocalVertsMap,
                            std::map<int,int>&&                    localV2globalV,
                            std::map<int,std::vector<int> >&&      Elem2Face,
                            std::map<int,std::vector<int> >&&      Elem2Vert,
                            std::map<int,std::vector<int> >&&      Face2Vert,
                            std::map<int,std::vector<int> >&&      Face2Elem,
                            std::map<int,int>                      partMap,
                            std::map<int,int>&&                    Elem2Rank,
                            std::set<int>                          TraceVertsOnRank,
                            std::set<int>                          TraceFacesOnRank,
                            std::set<int>&&                        ElemSet,
                            std::map<int,std::vector<int> >        ranges_id,
                            std::map<int,std::vector<int> >        ranges_ref);


            void buildUpdatedVertexAndFaceNumbering(MPI_Comm comm, 
                                                    std::map<int,int>               partMap,
                                                    std::set<int>                   TraceVertsOnRank,
                                                    std::set<int>                   TraceFacesOnRank,
                                                    std::map<int,std::vector<int> > ranges_id,
                                                    std::map<int,std::vector<int> > ranges_ref);



            void buildInteriorSharedAndBoundaryFaceMaps(MPI_Comm                    comm, 
                                                    std::map<int,int>               partMap,
                                                    std::set<int>                   TraceVertsOnRank,
                                                    std::set<int>                   TraceFacesOnRank,
                                                    std::map<int,std::vector<int> > ranges_id,
                                                    std::map<int,std::vector<int> > ranges_ref);
            
                int** getParMMGCommFace2GlobalVertMap();
                int** getParMMGCommFace2LocalVertMap();
                int* getParMMGCommColorFace();
                int* getParMMGCommNFacesPerColor();
                int getParMMGNComm();
                std::vector<int> getFace4ParMMG();
                std::map<int,int> getNewGlobalVert2LocalVertMap();
                std::map<int,int> getTagVert2GlobalVertMap();
                std::set<int> getElemSet();
                std::map<int,int> getGlobal2TagFMap();
                std::map<int,std::vector<double> > getLocalVertsMap();
                std::map<int,int> getLocalVert2GlobalVert();
                std::map<int,std::vector<int> > getElem2VertMap();
                std::map<int,std::vector<int> > getFace2VertMap();
                std::map<int,std::vector<int> > getFace2VertNewMap();
                std::map<int,int> getLocalSharedFace2GlobalSharedFace();

        private:
            // Data structures that are required for ParMMG.

                std::set<int>                          m_ElemSet;
                std::map<int,std::vector<int> >        m_Elem2Vert;
                std::map<int,std::vector<int> >        m_Elem2Face;
                std::map<int,std::vector<int> >        m_Face2Vert;
                std::map<int,std::vector<int> >        m_Face2Elem;
                std::map<int,int>                      m_locV2gloV;
                std::map<int,std::vector<double> >     m_LocalVertsMap;
                std::map<int,int>                      m_Elem2Rank;

                std::map<int,int>                       m_LocalV2GlobalV;
                std::map<int,int>                       m_SharedVertsOwned;
                std::map<int,int>                       m_NonSharedVertsOwned;
                std::map<int,int>                       m_SharedVertsNotOwned;
                std::map<int,int>                       m_sharedFace2RankMap;
                std::map<int,int>                       m_sharedVertexMapUpdatedGlobalID;
                std::map<int,int>                       m_sharedFaceMapUpdatedGlobalID;
                std::map<int,int>                       m_sharedVertex2RankMap;   
                std::map<int,int>                       m_boundaryFaces;
                std::map<int,int>                       m_boundaryFaces2Ref;
                std::map<int,std::vector<int> >         m_SharedFace2Rank;

                std::map<int,int>                       m_tag2globV;
                std::map<int,int>                       m_tag2globF;
                std::map<int,int>                       m_glob2tagF;
                std::map<int,int>                       m_glob2tagV;
                std::map<int,int>                       m_lhp;
                std::map<int,int>                       m_rhp;
                std::map<int,int>                       m_tagE2gE;
                std::map<int,int>                       m_gE2tagE;
                std::map<int,int>                       m_tagE2lE;
                std::map<int,int>                       m_lE2tagE;
                std::map<int,int>                       m_loc2globF;
                std::map<int,int>                       m_glob2locF;
                std::map<int,int>                       m_loc2globV;
                std::map<int,int>                       m_glob2locV;

                std::map<int,std::vector<int> >         m_Face2Elem_New;
                std::map<int,std::vector<int> >         m_Elem2Face_New;
                std::map<int,std::vector<int> >         m_Face2Vert_New;
                std::map<int, std::vector<int> >        m_sharedFace2Nodes;
                std::map<int, std::vector<int> >        m_colorRh;
                std::map<int, std::vector<int> >        m_colorFh;
                std::map<int, std::vector<int> >        m_interiorFace2Nodes;
                std::map<int, std::vector<int> >        m_boundaryFace2Nodes;
                std::map<int, std::vector<int> >        m_zone2bcface;
                
                std::map<int, std::vector<int> >        m_ColorsFaces;
                std::vector<int>                        m_faces4parmmg;
                std::map<int,int>                       m_globShF2locShF;
                std::map<int,int>                       m_locShF2globShF;

                int** m_ifc_tria_glo;
                int** m_ifc_tria_loc;
                int* m_color_face;
                int* m_ntifc;
                int m_ncomm;

                int Ne_glob;

};

#endif