#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate_lite.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"
#include "adapt_distri_parstate.h"
#include "adapt_prismtetratrace.h"

#ifndef ADAPT_REPARTITION_H
#define ADAPT_REPARTITION_H

class RepartitionObject{
        public:
                RepartitionObject(){};
                RepartitionObject(mesh* meshInput,
                                  std::map<int,std::vector<int> > elements2verts,
                                  std::map<int,std::vector<int> > elements2faces,
                                  std::map<int,std::vector<int> > elements2elements,
                                  PrismTetraTrace* trace,
                                  std::map<int,std::vector<double> > data,
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
                                                    std::map<int,std::vector<double> > data,
                                                    MPI_Comm comm, 
                                                    std::map<int,std::vector<int> >& elements2verts_update, 
                                                    std::map<int,std::vector<int> >& elements2faces_update, 
                                                    std::map<int,std::vector<int> >& elements2elements_update,
                                                    std::map<int,std::vector<double> >& elements2data_update);

                void DeterminePartitionLayout(std::map<int,std::vector<int> > elements, std::vector<int> element2rank, MPI_Comm comm);

                
                void DetermineElement2ProcMap(std::map<int,std::vector<int> >     ien, 
                                                 std::map<int,std::vector<int> >   ief,
                                                 std::map<int,std::vector<int> >   iee,  
                                                 std::map<int,std::vector<double> > data,
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

                
                std::map<int,std::vector<double> > getElement2DataMap();

                std::map<int, std::vector<int> > getGlobalElement2LocalVertMap();

                std::vector<std::vector<double> > getLocalVerts();

                std::map<int, std::vector<double> > getLocalVertsMap();

                std::map<int,int> getGlobalElement2Rank();

                std::map<int,std::vector<int> > getElement2VertexMap();

                std::map<int,std::vector<int> > getElement2ElementMap();

                std::map<int,std::vector<int> > getElement2FacesMap();

                std::map<int, std::vector<int> > getFace2VertexMap();

                std::map<int, std::vector<int> > getFace2NVertexMap();

                std::vector<int> getLocElem();
                        
                void getAdjacentElementLayer(std::map<int,std::vector<int> > element2verts,
                                             std::map<int,std::vector<int> > element2faces,
                                             std::map<int,std::vector<int> > element2element,
                                             PrismTetraTrace* trace,
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

        private:
                std::vector<int> part;
                std::map<int,int> part_global;
                std::map<int,int> part_map;

                int eloc;
                int vloc;
                int floc;


                std::map<int, std::vector<int> > globElem2globVerts;
                std::map<int, std::vector<int> > globElem2locVerts;

                std::vector<std::vector<int> > lelement2lvertex;

                std::map<int,int> lvertex2gvertex;
                std::map<int,int> gvertex2lvertex;

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
                //std::vector<std::vector<double> > LocalVerts;
                std::map<int, std::vector<double> > LocalVertsMap;
                int itel_locadj;
                int nLoc_Elem;
                int nLocAndAdj_Elem;
                int nLoc_Verts;
                std::map<int,std::vector<int> > globElem2localFaces;
                std::map<int,std::vector<int> > globElem2globFaces;
                std::map<int,std::vector<int> > globFace2GlobalElements;
                std::vector<int> Loc_Elem;
                std::vector<int> Loc_Elem_Nv;
                std::vector<int> Loc_Elem_Nf;

                std::map<int,int> LocElem2Nv;
                std::map<int,int> LocElem2Nf;
                std::map<int,int> LocElem2Ndata;
                std::vector<int> LocAndAdj_Elem;
                std::vector<int> LocAndAdj_Elem_Nv;
                std::vector<int> LocAndAdj_Elem_Nf;
                std::map<int,int> LocalElement2GlobalElement;
                std::map<int,int> GlobalElement2LocalElement;
                std::map<int, std::vector<int> > globVerts2globElem;

                std::map<int, std::vector<int> > elements2verts_update;
                std::map<int, std::vector<int> > elements2faces_update;
                std::map<int, std::vector<int> > elements2elements_update;
                std::map<int, std::vector<double> > elements2data_update;
                std::map<int, std::vector<int> > face2verts_update;
                std::map<int, std::vector<int> > face2Nverts_update;

};      

#endif


