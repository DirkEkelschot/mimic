#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"
#include "adapt_distri_parstate.h"

#ifndef ADAPT_PARTITION_H
#define ADAPT_PARTITION_H

class Partition {
   public:
    Partition(){};
    Partition(ParArray<int>* ien, ParArray<int>* iee, ParArray<int>* ief, ParArray<int>* ie_Nv, ParArray<int>* ieie_Nf, ParArray<int>* ifn, ParArray<int>* ife, ParArray<int>* if_ref, ParArray<int>* if_Nv,  ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, ParallelState* ife_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, Array<int>* tetCnt, MPI_Comm comm);
    ~Partition();
    void DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, MPI_Comm comm);
    void DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* ie_Nv, ParArray<int>* ie_Nf, ParArray<double>* xcn, Array<double>* U, Array<int>* tetCnt, MPI_Comm comm);
    void DetermineAdjacentElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* part, ParArray<double>* xcn, Array<double>* U, MPI_Comm comm);
    
    std::vector<int> getAdjacentElementLayer(ParArray<int>* ien,
                                 std::vector<int> Loc_Elem_input,
                                 std::vector<int> Loc_Elem_Nf_input,
                                            std::map<int,std::vector<int> > iee_vec, ParArray<int>* part, ParArray<double>* xcn, Array<double>* U, MPI_Comm comm);
    
    std::vector<int> UpdateAdjacentElementLayer(ParArray<int>* ien,
                                                std::vector<int> Loc_Elem_Packed,
                                            std::map<int,std::vector<int> > iee_vec, ParArray<int>* part, ParArray<double>* xcn, Array<double>* U, MPI_Comm comm);
    
    void CreatePartitionDomain();
    void CreatePartitionDomainTest();
    std::map<int,Array<double>* > PartitionAuxilaryData(Array<double>* U, MPI_Comm comm);
    std::map<int,double> CommunicateLocalDataUS3D(Array<double>* U, MPI_Comm comm);
    void AddStateForAdjacentElements(std::map<int,double> U, MPI_Comm comm);
    void AddStateVecForAdjacentElements(std::map<int,Array<double>* > &U, int nvar, MPI_Comm comm);
    void AddAdjacentVertexDataUS3D(std::map<int,double> &Uv, MPI_Comm comm);
    void AddStateVecForAdjacentVertices(std::map<int,Array<double>* > &Uv, int nvar, MPI_Comm comm);
    i_part_map* getElement2EntityPerPartition(ParArray<int>* iee, std::vector<int> Loc_Elem,  std::vector<int> Loc_Elem_Ne, MPI_Comm comm);
    i_part_map* UpdateElement2EntityPerPartition(ParArray<int>* iee, std::vector<int> LocAndAdj_Elem, std::vector<int> LocAndAdj_Elem_Ne, MPI_Comm comm);
    void UpdateElement2EntityPerPartition_V2(ParArray<int>* iee, std::vector<int> LocAndAdj_Elem_Packed, int index, i_part_map* iee_part_map_input, MPI_Comm comm);
    void UpdateFace2EntityPerPartition_V2(ParArray<int>* ife, std::vector<int> LocAndAdj_Elem_Packed, int index, i_part_map* ief_part_map_input, i_part_map* ifn_part_map_input, MPI_Comm comm);
    i_part_map* getFace2EntityPerPartition(i_part_map* ief_part_map_input, ParArray<int>* ife, MPI_Comm comm);
    i_part_map* getFace2NodePerPartition(ParArray<int>* ifn, MPI_Comm comm);
    Domain* getPartitionDomain();
    std::map<int,double> ReduceFieldToVertices(std::map<int,double> Uelem);
    std::map<int,double> ReduceFieldToAllVertices(std::map<int,double> Uelem);
    std::map<int,Array<double>* > ReduceStateVecToAllVertices(std::map<int,Array<double>* > UaddAdj, int nvar);
    std::map<int,Array<double>* > ReduceStateVecToAllVertices_V2(std::map<int,Array<double>* > UaddAdj, int nvar);
    std::map<int,Array<double>*> ReduceMetricToVertices(std::map<int,Array<double>* > Telem);
    std::map<int,int> getGlobalVert2GlobalElement();
    void GetPartitionFaceColorMaps();
    std::vector<int> getTetCnt();
    std::vector<int> getLocElem();
    std::vector<int> getLocElemNv();
    std::map<int,int> getLocElem2Nv();
    std::map<int,int> getLocElem2Nf();
    std::vector<double> getLocElemVaria();
    std::map<int,Array<double>* > getLocAndAdjElemVaria();
    int getnLoc_Elem();
    std::vector<int> getLocAndAdjElem();
    std::vector<int> getLocAndAdjElem_Nv();
    std::vector<int> getLocAndAdjElem_Nf();
    int getnLocAndAdj_Elem();
    int getNloc_Elem();
    int getNLocAndAdj_Elem();
    int getnLoc_Verts();
    int* getXadj();
    int* getAdjcny();
    std::vector<int> getLocAndAdj_Elem();
    std::vector<int> getLocAndAdj_Elem_Nv();
    ParArray<int>* getLocalPartition();
    Array<int>* getGlobalPartition();
    std::vector<Vert*> getLocalVerts();
    std::map<int,std::map<int,double> > getNode2NodeMap();
    std::map<int,std::map<int,double> > getNode2NodeMap_V2();
    std::map<int,std::set<int> > getSecondLayerAdjacency();
    Vert* getLocalVert(int v_loc_id);
    
    std::vector<std::vector<int> > getLocalElem2GlobalVert();
    std::vector<std::vector<int> > getLocalElem2LocalVert();
    
    std::map<int,int> getLocalVert2GlobalVert();
    std::map<int,int> getGlobalVert2LocalVert();

    
    std::map<int,int> getLocalElement2GlobalElement();
    std::map<int,int> getGlobalElement2LocalElement();

    std::map<int,int> getLocalFace2GlobalFace();
    std::map<int,int> getGlobalFace2LocalFace();
    std::map<int,std::vector<int> > getglobElem2localFaces();
    std::map<int,std::vector<int> > getglobElem2globFaces();
    std::map<int,std::vector<int> > getglobFace2GlobalElements();
    std::map<int,std::vector<int> > getGlobElem2GlobVerts();
    std::map<int,std::vector<int> > getGlobElem2LocVerts();
    std::map<int,std::vector<int> > getGlobVert2GlobElem();
    std::set<int> getElemSet();
    std::set<int> getLocElemSet();
    std::vector<double> getUelem();
    double getU0atGlobalElem(int elem);
    double getUauxatGlobalElem(int elem);
    Array<double>* getUvert();
    ParallelState* getXcnParallelState();
    ParallelState* getIenParallelState();
    ParallelState* getParallelState();   
    ParallelState_Parmetis* getParallelStateParmetis();
    
    i_part_map* getIEEADJ2partmap();
    i_part_map* getIEEADJpartmap();
    i_part_map* getIEEpartmap();
    
    i_part_map* getIEFADJ2partmap();
    i_part_map* getIEFADJpartmap();
    i_part_map* getIEFpartmap();

    i_part_map* getIFNADJ2partmap();
    i_part_map* getIFNADJpartmap();
    i_part_map* getIFNpartmap();
    
    i_part_map* getIENpartmap();

    i_part_map* getIFERankpartmap();
    i_part_map* getIF_Nvpartmap();
    i_part_map* getIFEpartmap();
    i_part_map* getIFREFpartmap();
    i_part_map* getIE_Nfpartmap();
    
    std::map<int,std::map<int,double> > getNode2Element(i_part_map* iee_part_map_input, std::map<int,Vert*> ElemCentroids, std::map<int,double> ElemVolumes);
    std::map<int,std::map<int,double> > getNode2Element_V2(i_part_map* iee_part_map_input, std::map<int,Vert*> ElemCentroids, std::map<int,double> ElemVolumes);
    
   private:
      MPI_Comm comm_p;
      int rank;
      int size;
      std::vector<int> Loc_Elem;
      std::vector<int> Loc_Elem_Nv;
      std::vector<int> Loc_Elem_Nf;
      std::vector<double> Loc_Elem_Varia;
      std::vector<int> Loc_Elem_TetCnt;
      std::vector<int> LocAndAdj_Elem;
      std::vector<int> LocAndAdj_Elem_Nv;
      std::vector<int> LocAndAdj_Elem_Nf;
      std::vector<int> LocAndAdj_Elem_Varia;
      std::map<int,int> LocElem2Nv;
      std::map<int,int> LocElem2Nf;
      std::map<int,Array<double>* > LocAndAdjElemVaria;
      int NelGlob;
      int nloc;
      int eloc;
      int vloc;
      int floc;
      int* xadj;
      int* adjcny;
      
      std::vector<int> loc_r_elem;
      std::vector<int> loc_r_nv_elem;
      std::vector<int> loc_r_nf_elem;
      std::vector<double> loc_varia;
      std::vector<int> loc_tetc;
      Domain* pDom;
      int itel_locadj;
      int nLoc_Elem;
      int nLocAndAdj_Elem;
      int nLoc_Verts;
      std::set<int> elem_set;
      std::map<int,int> elem_map;
      std::set<int> loc_r_elem_set;
      //Array<int>* LocAndAdj_Elem;
      ParArray<int>* part;
      Array<int>* part_global;
      std::vector<Vert*> LocalVerts;

      std::set<int> unique_vertIDs_on_rank_set;
      std::vector<int> unique_verts_on_rank_vec;
      std::set<int> unique_faceIDs_on_rank_set;
    
      std::map<int, std::vector<int> > globVerts2globElem;
    
      std::map<int, std::vector<int> > globElem2globVerts;
      std::map<int, std::vector<int> > globElem2locVerts;
      std::vector<std::vector<int> > LocalElem2GlobalVert;
      std::vector<std::vector<int> > LocalElem2LocalVert;
      std::map<int,int> LocalVert2GlobalVert;
      std::map<int,int> GlobalVert2LocalVert;
    

      std::map<int,int> LocalFace2GlobalFace;
      std::map<int,int> GlobalFace2LocalFace;
      std::map<int,std::vector<int> > globElem2localFaces;
      std::map<int,std::vector<int> > globElem2globFaces;
      std::map<int,std::vector<int> > globFace2GlobalElements;

      std::map<int,int> LocalElement2GlobalElement;
      std::map<int,int> GlobalElement2LocalElement;
    
      //Array<double>* U0Elem; // This is the value of U0 for each cell.
      Array<double>* U0Vert; // This is the reduced average for each vert based on
      std::map<int,std::vector<double> > collect_var;
      ParallelState* xcn_pstate;
      ParallelState* ien_pstate;
      ParallelState* ife_pstate;
      ParallelState_Parmetis* pstate_parmetis;
    
    
      std::map<int,std::vector<int> > adj_elements;
      ScheduleObj* adj_schedule;
      ScheduleObj* part_schedule;
      std::map<int,std::vector<int> > elms_to_send_to_ranks;
      std::map<int,std::vector<int> > nvPerElms_to_send_to_ranks;
      std::map<int,std::vector<int> > nfPerElms_to_send_to_ranks;
      std::map<int,std::vector<int> > part_tot_recv_elIDs;
      std::map<int,std::vector<double> > part_tot_recv_varias;
      std::map<int,std::vector<int> > reqstd_adj_ids_per_rank;
      std::map<int,std::map<int,Vert*> > node2elem_map;
      std::map<int,std::map<int,double> > node2elemVol_map;

      i_part_map* if_Erank_part_map;
      i_part_map* if_Nv_part_map;
      i_part_map* ifn2_part_map;
      i_part_map* if_ref_part_map;
      i_part_map* ifn_part_map;
      i_part_map* ifn_adj_part_map;
      i_part_map* ifn_adj2_part_map;
      i_part_map* ife_part_map;
      i_part_map* ieNf_part_map;
      i_part_map* iee_part_map;
      i_part_map* iee_adj_part_map;
      i_part_map* iee_adj2_part_map;
      i_part_map* ief_part_map;
      i_part_map* ief_adj_part_map;
      i_part_map* ief_adj2_part_map;
      i_part_map* ien_part_map;
      i_part_map* ien_adj_part_map;
      i_part_map* ien_adj2_part_map;
};
#endif
