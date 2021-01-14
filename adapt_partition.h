#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"

#ifndef ADAPT_PARTITION_H
#define ADAPT_PARTITION_H

class Partition {
   public:
    Partition(){};
    Partition(ParArray<int>* ien, ParArray<int>* iee, ParArray<int>* ief, ParArray<int>* ifn, ParArray<int>* ife, ParArray<int>* if_ref,  ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, ParallelState* ifn_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm);
    
    void DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, MPI_Comm comm);
    void DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParallelState* ien_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm);
    void DetermineAdjacentElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* part, ParallelState* ien_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm);
    void DetermineAdjacentElement2ProcMapUS3D(ParArray<int>* ien, std::map<int,std::vector<int> > iee_vec, ParArray<int>* part, ParallelState* ien_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm);
    std::vector<double> PartitionAuxilaryData(Array<double>* U, MPI_Comm comm);
    std::map<int,double> CommunicateLocalDataUS3D(Array<double>* U, MPI_Comm comm);
    std::map<int,double> CommunicateAdjacentDataUS3D(Array<double>* U, MPI_Comm comm);
    
    i_part_map* getElement2EntityPerPartition(ParArray<int>* iee, ParallelState* ien_pstate, MPI_Comm comm);
    i_part_map* getFace2EntityPerPartition(ParArray<int>* ife, ParallelState* ife_pstate, MPI_Comm comm);
    std::vector<int> getLocElem();
    int getnLoc_Elem();
    std::vector<int> getLocAndAdjElem();
    int getnLocAndAdj_Elem();
    int getNloc_Elem();
    int getNLocAndAdj_Elem();
    int getnLoc_Verts();
    int* getXadj();
    int* getAdjcny();
    std::vector<int> getLocAndAdj_Elem();
    ParArray<int>* getLocalPartition();
    Array<int>* getGlobalPartition();
    std::vector<Vert> getLocalVerts();
    Vert getLocalVert(int v_loc_id);
    
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
    
    i_part_map* getIEEpartmap();
    i_part_map* getIEFpartmap();
    i_part_map* getIENpartmap();
    
    i_part_map* getIFNpartmap();
    i_part_map* getIFEpartmap();
    i_part_map* getIFREFpartmap();
    
   private:
      
      std::vector<int> Loc_Elem;
      std::vector<int> LocAndAdj_Elem;
      int nloc;
      int eloc;
      int vloc;
      int floc;
      int* xadj;
      int* adjcny;
      std::vector<int> loc_elem;
      std::vector<double> loc_rho;
      int nLoc_Elem;
      int nLocAndAdj_Elem;
      int nLoc_Verts;
      std::set<int> elem_set;
      std::set<int> loc_elem_set;
      //Array<int>* LocAndAdj_Elem;
      ParArray<int>* part;
      Array<int>* part_global;
      std::vector<Vert> LocalVerts;
      

      std::set<int> unique_vertIDs_on_rank_set;
      std::vector<int> unique_verts_on_rank_vec;
      std::set<int> unique_faceIDs_on_rank_set;
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
      ParallelState_Parmetis* pstate_parmetis;
    
    
      std::map<int,std::vector<int> > adj_elements;
      ScheduleObj* adj_schedule;
      ScheduleObj* part_schedule;
      std::map<int,std::vector<int> > elms_to_send_to_ranks;
      std::map<int,std::vector<int> > part_tot_recv_elIDs;
      std::map<int,std::vector<int> >  reqstd_adj_ids_per_rank;
      
      i_part_map* if_ref_part_map;
      i_part_map* ifn_part_map;
      i_part_map* ife_part_map;
    
      i_part_map* iee_part_map;
      i_part_map* ief_part_map;
      i_part_map* ien_part_map;
};
#endif
