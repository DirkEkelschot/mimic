#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate.h"
#include "adapt_operations.h"

#ifndef ADAPT_PARTITION_H
#define ADAPT_PARTITION_H

class Partition {
   public:
    Partition(ParArray<int>* ien, ParArray<int>* ief, ParallelState_Parmetis* pstate_parmetis, ParallelState* pstate, ParArray<double>* xcn, ParallelState* xcn_parstate, ParArray<double>* U, MPI_Comm comm);
    void DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, ParallelState* pstate, MPI_Comm comm);
    void DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* part, ParallelState* pstate, ParArray<double>* xcn, ParallelState* xcn_parstate, ParArray<double>* U, MPI_Comm comm);
    int getNlocElem();
    int getNlocVerts();
    int* getXadj();
    int* getAdjcny();
    ParArray<int>* getPart();
    std::vector<Vert> getLocalVerts();
    Vert getLocalVert(int v_loc_id);
    
    Array<int>* getLocalElem2GlobalVert();
    Array<int>* getLocalElem2LocalVert();
    
    std::map<int,int> getLocalVert2GlobalVert();
    std::map<int,int> getGlobalVert2LocalVert();
    
    Array<int>* getLocalElem2GlobalFace();
    Array<int>* getLocalElem2LocalFace();
    
    std::map<int,int> getLocalElement2GlobalElement();
    std::map<int,int> getGlobalElement2LocalElement();

    std::map<int,int> getLocalFace2GlobalFace();
    std::map<int,int> getGlobalFace2LocalFace();
    std::map<int,std::vector<int> > getElem2GlobalFace();
    std::map<int,std::vector<int> > getGlobalFace2Elem();
    std::map<int,std::vector<int> > getGlobElem2GlobVerts();
    std::map<int,std::vector<int> > getGlobElem2LocVerts();
    std::set<int> getElemSet();
    Array<double>* getUelem();
    double getU0atGlobalElem(int elem);
    Array<double>* getUvert();
    
    ParallelState* getXcnParallelState();
    ParallelState* getIenParallelState();
    
    ParallelState_Parmetis* getParallelStateParmetis();
    
   private:
      
      int* xadj;
      int* adjcny;
    
      int NlocElem;
      int NlocVerts;
      std::set<int> elem_set;
      Array<int>* ElemPart;
      ParArray<int>* part;
      Array<int>* part_global;
      std::vector<Vert> LocalVerts;
    
        
      std::map<int, std::vector<int> > globElem2globVerts;
      std::map<int, std::vector<int> > globElem2locVerts;
      Array<int>* LocalElem2GlobalVert;
      Array<int>* LocalElem2LocalVert;
      std::map<int,int> LocalVert2GlobalVert;
      std::map<int,int> GlobalVert2LocalVert;
    
      Array<int>* LocalElem2GlobalFace;
      Array<int>* LocalElem2LocalFace;
      std::map<int,int> LocalFace2GlobalFace;
      std::map<int,int> GlobalFace2LocalFace;
      std::map<int,std::vector<int> > Elem2GlobalFace;
      std::map<int,std::vector<int> > GlobalFace2Elem;

      std::map<int,int> LocalElement2GlobalElement;
      std::map<int,int> GlobalElement2LocalElement;
    
      Array<double>* U0Elem; // This is the value of U0 for each cell.
      Array<double>* U0Vert; // This is the reduced average for each vert based on
    
      ParallelState* xcn_pstate;
      ParallelState_Parmetis* pstate_parmetis;
};


inline Partition::Partition(ParArray<int>* ien, ParArray<int>* ief, ParallelState_Parmetis* pstate_parmetis, ParallelState* pstate, ParArray<double>* xcn, ParallelState* xcn_parstate, ParArray<double>* U, MPI_Comm comm)
{
    
    // This function computes the xadj and adjcny array and the part array which determines which element at current rank should be sent to other ranks.
    DeterminePartitionLayout(ien, pstate_parmetis, pstate, comm);
    
    // This function takes care of the send and receive operations in order to send the appropriate elements and corresponding vertices to the appropriate rank.
    // These operations are based on the fact that part holds the desired partitioning of the elements. the spread of the vertices is based on the fact that all the vertices stored in xcn are distributed "uniformly";
    DetermineElement2ProcMap(ien, ief, part, pstate, xcn, xcn_parstate, U, comm);
}

inline void Partition::DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, ParallelState* pstate, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = ien->getNrow();
    int nloc = nrow;

    //=================================================================
    //=================================================================
    //=================================================================
    
    //ParallelState_Parmetis* pstate_parmetis2 = new ParallelState_Parmetis(ien,comm,8);
//
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.1};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt = NULL;

    int np           = size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    int* part_arr = new int[nloc];
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);
    
    
    for(int u=0;u<nloc;u++)
    {
        part_arr[u] = rank;
    }
    /*
    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(),
                           xadj_par, adjncy_par,
                               elmwgt, adjwgt,
                       vsize, wgtflag,
                   numflag, ncon, nparts,
                   tpwgts, ubvec, itr, options,
                   &edgecut, part_arr, &comm);
    */
    /*
    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);
    
    */
    part = new ParArray<int>(ien->getNglob(),1,comm);
    part_global = new Array<int>(ien->getNglob(),1);

    part->data = part_arr;
    xadj = xadj_par;
    adjcny = adjncy_par;
    
    MPI_Allgatherv(&part->data[0],
                   nloc, MPI_INT,
                   &part_global->data[0],
                   pstate->getNlocs(),
                   pstate->getOffsets(),
                   MPI_INT,comm);
    
}


inline void Partition::DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* part, ParallelState* pstate, ParArray<double>* xcn, ParallelState* xcn_parstate, ParArray<double>* U, MPI_Comm comm)
{
    
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int el_id;
    int p_id;
    int v_id;
    Vert V;
    std::vector<Vert> part_verts;
    std::vector<std::vector<int> > part_elem2verts;

    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceIDs_to_send_to_ranks;
    std::map<int,std::vector<double> > rho_to_send_to_ranks;

    std::set<int> unique_vertIDs_on_rank_set;
    std::vector<int> unique_verts_on_rank_vec;
    std::set<int> unique_faceIDs_on_rank_set;
    
    std::set<int> u_verts_part_set;
    std::vector<int> u_verts_part_vec;
    std::map<int,std::vector<int> > rank2req_vert;

    std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::vector<int> part_v;
    int r     = 0;
    int lv_id = 0;
    int f_id  = 0;
    std::map<int,int> part_uvert_loc2glob;
    std::map<int,int> part_uvert_glob2loc;
    
    std::map<int, int> part_fuvert_loc2glob;
    std::map<int, int> part_fuvert_glob2loc;
    int xcn_o = xcn->getOffset(rank);
    int xcn_n = xcn->getNloc(rank);
    int vloc = 0;
    int floc = 0;
    std::vector<int> loc_elem;
    std::vector<double> loc_rho;
    int ien_o = part->getOffset(rank);
    double rho = 0.0;
    int not_on_rank=0;
    int on_rank = 0;
    std::set<int> requested_elements;
    int* new_offsets = new int[size];
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = xcn_parstate->getOffsets()[i]-1;
    }
            
    for(i=0;i<part->getNrow();i++)
    {
        p_id  = part->getVal(i,0);
        el_id = part->getOffset(rank)+i;
        rho   = U->getVal(i,0);
        //std::cout << "rho = " << rho << std::endl;
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            rho_to_send_to_ranks[p_id].push_back(rho);
            
            //====================Hexes=======================
            for(int k=0;k<8;k++)//This works for hexes.
            {
                v_id = ien->getVal(i,k);
                vertIDs_to_send_to_ranks[p_id].push_back(v_id);
                
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end())// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_offsets,size,v_id);

                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vloc++;
                    }
                }
                
                if (k<6)// faces
                {
                    f_id = ien->getVal(i,k);
                    faceIDs_to_send_to_ranks[p_id].push_back(f_id);
                }
                
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            //====================Hexes=======================
            
            not_on_rank++;
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;
            for(int k=0;k<8;k++)// looping over the vertices for element "i".
            {
                v_id = ien->getVal(i,k);
                
                //elem.push_back(v_id);
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end())// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_offsets,size,v_id);

                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vloc++;
                    }
                    lv_id++;
                }
                
                if(k<6) // just store all faceID since they map to the local elemID.
                {
                    f_id = ief->getVal(i,k);
                    if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
                    {
                        unique_faceIDs_on_rank_set.insert(f_id);
                        faceIDs_on_rank.push_back(f_id);
                        floc++;
                    }
                }
            }
            on_rank++;
        }
        
        loc_elem.push_back(el_id);
        loc_rho.push_back(rho);
        elem_set.insert(el_id);
    }
    
    int nRank_reqElems = elms_to_send_to_ranks.size();
    
    int* reduced_nRank_reqElems = new int[size];
    int* arr_nRank_reqElems = new int[size];
    
    for(i=0;i<size;i++)
    {
        reduced_nRank_reqElems[i] = 0;
        
        if(i==rank)
        {
            arr_nRank_reqElems[i] = nRank_reqElems+1;
        }
        else
        {
            arr_nRank_reqElems[i] = 0;
        }
    }

    MPI_Allreduce(arr_nRank_reqElems, reduced_nRank_reqElems, size, MPI_INT, MPI_SUM, comm);
    
    int* reduced_nRank_reqElems_offset = new int[size];
    int offset = 0;
    for(i=0;i<size;i++)
    {
        reduced_nRank_reqElems_offset[i] = offset;
        offset = offset+reduced_nRank_reqElems[i];
    }
    int nTot_reqElements = 0;
    int nRank_reqElems_p_one = nRank_reqElems+1; // This size is added by one since we add the rank number to the array.
    
    MPI_Allreduce(&nRank_reqElems_p_one, &nTot_reqElements, 1, MPI_INT, MPI_SUM, comm);

    int* sendFromRank2Rank_elem_Global    = new int[nTot_reqElements];
    int* sendNelemFromRank2Rank_Global    = new int[nTot_reqElements];
    
    for(i=0;i<nTot_reqElements;i++)
    {
        sendFromRank2Rank_elem_Global[i]  = 0;
        sendNelemFromRank2Rank_Global[i]  = 0;
    }
    
    
    int* send2Rank_fromRank_e       = new int[nRank_reqElems+1];
    int* sendNelem2Rank_from_rank   = new int[nRank_reqElems+1];
    
    send2Rank_fromRank_e[0]         = rank;
    sendNelem2Rank_from_rank[0]     = -1;
    int t = 1;
    
    std::map<int,std::vector<int> >::iterator it;
    for(it=elms_to_send_to_ranks.begin();it!=elms_to_send_to_ranks.end();it++)
    {
        send2Rank_fromRank_e[t]          = it->first;
        sendNelem2Rank_from_rank[t]      = it->second.size();
        t++;
    }

    MPI_Allgatherv(&send2Rank_fromRank_e[0],
                   nRank_reqElems_p_one, MPI_INT,
                   &sendFromRank2Rank_elem_Global[0],
                   reduced_nRank_reqElems,
                   reduced_nRank_reqElems_offset,
                   MPI_INT,comm);
    
    MPI_Allgatherv(&sendNelem2Rank_from_rank[0],
                   nRank_reqElems_p_one, MPI_INT,
                   &sendNelemFromRank2Rank_Global[0],
                   reduced_nRank_reqElems,
                   reduced_nRank_reqElems_offset,
                   MPI_INT,comm);
    
    
    std::map<int, std::set<int> > sendFromRank2Rank_e;
    std::map<int, std::vector<int> > recvRankFromRank_map_e;
    std::map<int, std::vector<int> > recvRankFromRank_map_Nelem_e;
    
    //=========================================
    for(i=0;i<size;i++)
    {
        int of = reduced_nRank_reqElems_offset[i];
        int nl = reduced_nRank_reqElems[i];
        for(int j=of+1;j<of+nl;j++)
        {
            sendFromRank2Rank_e[sendFromRank2Rank_elem_Global[of]].insert(sendFromRank2Rank_elem_Global[j]);
            recvRankFromRank_map_e[sendFromRank2Rank_elem_Global[j]].push_back(sendFromRank2Rank_elem_Global[of]);
            recvRankFromRank_map_Nelem_e[sendFromRank2Rank_elem_Global[j]].push_back(sendNelemFromRank2Rank_Global[j]);
        }
    }

    
    std::map<int,std::map<int,int> > RecvAlloc_map_e;
    for(it=recvRankFromRank_map_e.begin();it!=recvRankFromRank_map_e.end();it++)
    {
        for(int k=0;k<it->second.size();k++)
        {
            RecvAlloc_map_e[it->first][recvRankFromRank_map_e[it->first][k]] = recvRankFromRank_map_Nelem_e[it->first][k];
        }
    }
    
    //============================================================================
    //================ Print the receiving schedule for now;======================
    //============================================================================
//    if(rank == 0)
//    {
//        for(it=recvRankFromRank_map_e.begin();it!=recvRankFromRank_map_e.end();it++)
//        {
//            std::cout << "rank " << it->first << " receives an array of size ";
//            for(int k=0;k<it->second.size();k++)
//            {
//                std::cout << recvRankFromRank_map_Nelem_e[it->first][k] << " from " << recvRankFromRank_map_e[it->first][k] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    //============================================================================
    //================ Print the receiving schedule for now;======================
    //============================================================================
    
    
    std::map<int,int> alloc = RecvAlloc_map_e[rank];
    
    int* recv_offset = new int[RecvAlloc_map_e.size()+1];
    recv_offset[0]   = 0;
    int* recv_loc    = new int[RecvAlloc_map_e.size()];;
    int TotNelem_recv    = 0;
    
    std::map< int, int> RecvAlloc_offset_map_e;
    std::map< int, int>::iterator it_loc;
    i = 0;
    for(it_loc=alloc.begin();it_loc!=alloc.end();it_loc++)
    {
        recv_loc[i]      = it_loc->second;
        recv_offset[i+1] = recv_offset[i]+recv_loc[i];
        TotNelem_recv    = TotNelem_recv+it_loc->second;
        RecvAlloc_offset_map_e[it_loc->first]=recv_offset[i];

        i++;
    }
    
    int* TotRecvElement_IDs       = new int[TotNelem_recv];
    double* TotRecvElement_rhos   = new double[TotNelem_recv];

    int* TotRecvElement_IDs_v     = new int[TotNelem_recv*8];
    int* TotRecvElement_IDs_f     = new int[TotNelem_recv*6];
    
    for(int i =0;i<TotNelem_recv;i++)
    {
        TotRecvElement_IDs[i] = 0;
        TotRecvElement_rhos[i] = 0;
    }
    
    int n_req_recv;
    
    int n_req_recv_v;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = elms_to_send_to_ranks.begin(); it != elms_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = n_req*8;
                int n_req_f         = n_req*6;
                int dest            = it->first;
                                
                MPI_Send(&n_req, 1, MPI_INT, dest, dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                MPI_Send(&rho_to_send_to_ranks[it->first][0], n_req, MPI_DOUBLE, dest, 20000+100+dest*2, comm);
                
                //MPI_Send(&n_req_v, 1, MPI_INT, dest, 9000+dest, comm);
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                i++;
            }
        }
        else if (sendFromRank2Rank_e[q].find( rank ) != sendFromRank2Rank_e[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_IDs[RecvAlloc_offset_map_e[q]], n_req_recv, MPI_INT, q, 100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_rhos[RecvAlloc_offset_map_e[q]], n_req_recv, MPI_DOUBLE, q, 20000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            //MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, 9000+rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_IDs_v[RecvAlloc_offset_map_e[q]*8], n_req_recv*8, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_IDs_f[RecvAlloc_offset_map_e[q]*6], n_req_recv*6, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
        }
    }
    
    int Nel_extra = TotNelem_recv;
    int cnt_v = 0;
    int cnt_f = 0;
    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;
        
        for(int k=0;k<8;k++)
        {
            int v_id_n = TotRecvElement_IDs_v[cnt_v];
            //elem.push_back(v_id_n);
            r = FindRank(new_offsets,size,v_id_n);
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id_n);
                //part_v.push_back(r);
                
                if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
                {
                    rank2req_vert[r].push_back(v_id_n); // add vertex to rank2req_vert map.
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n); // add the vertex to list that is already available on rank.
                    vloc++;
                }
            }
            cnt_v++;
            
            if(k<6)
            {
                int f_id_n = TotRecvElement_IDs_f[cnt_f];
                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
                {
                    unique_faceIDs_on_rank_set.insert(f_id);
                    faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                    floc++;
                }
                cnt_f++;

            }
            
        }
        //part_elem2verts.push_back(elem);
        elem_set.insert(TotRecvElement_IDs[i]);
    }
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    std::vector<int> adj_elements;
    for(int i=0;i<part->getNrow();i++)
    {
        int start = xadj[i];
        int end   = xadj[i+1];
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            p_id = part_global->getVal(adjEl_id,0);
            if(p_id!=rank && (elem_set.find(adjEl_id)==elem_set.end()))
            {
		adj_elements.push_back(adjEl_id);
		elem_set.insert(adjEl_id);
                req_elem[p_id].push_back(adjEl_id);
                itel++;
            }
        }
    }
    
    ScheduleObj* sobj_el = DoScheduling(req_elem,comm);
    
    std::map<int,std::vector<int> >  reqstd_adj_ids_per_rank;
    int n_reqstd_adj_ids;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = req_elem.begin(); it != req_elem.end(); it++)
            {
                int n_req_adj_el           = it->second.size();
                int dest                   = it->first;
                
                int destination = dest;
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req_adj_el, 1, MPI_INT, dest, 9876000+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req_adj_el, MPI_INT, dest, 9876000*2+dest*2, comm);
               i++;
            }
        }
        else if (sobj_el->SendFromRank2Rank[q].find( rank ) != sobj_el->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_adj_ids, 1, MPI_INT, q, 9876000+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recv_reqstd_adj_ids(n_reqstd_adj_ids);
            MPI_Recv(&recv_reqstd_adj_ids[0], n_reqstd_adj_ids, MPI_INT, q, 9876000*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_adj_ids_per_rank[q] = recv_reqstd_adj_ids;
                 
        }
    }
    
    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,std::vector<int> > send_adj_verts_IDs;
    std::map<int,std::vector<int> > send_adj_faces_IDs;
    int TotNelem_adj_recv = 0;
    std::vector<int> TotAdj_El_IDs;
    std::map<int,std::vector<int> > send_adj_rhos;
    std::vector<double> TotAdj_Rhos;
    int adj_id;

    int offset_new = pstate->getOffset(rank);
    int nloc_new   = pstate->getNloc(rank);

    
    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            adj_id = itv->second[j];
            TotAdj_El_IDs.push_back(adj_id);
            TotAdj_Rhos.push_back(U->getVal(adj_id-offset_new,0));
            send_adj_rhos[dest].push_back(U->getVal(adj_id-offset_new,0));

            for(int k=0;k<8;k++)
            {
                v_id = ien->getVal(adj_id-offset_new,k);
                send_adj_verts_IDs[dest].push_back(v_id);
            }
            
            for(int k=0;k<6;k++)
            {
                int offset_new = pstate->getOffset(rank);
                f_id = ief->getVal(adj_id-offset_new,k);
                send_adj_faces_IDs[dest].push_back(f_id);
            }
        }
        
        TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    //std::cout << "itel = " << itel << " " << TotNelem_adj_recv << " " << TotAdj_El_IDs.size() << std::endl;
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn = 0;
    std::map<int,int  > recv_adj_back_Nverts;
    std::map<int,int* > recv_adj_back_verts_ids;
    std::map<int,int  > recv_adj_back_Nfaces;
    std::map<int,int* > recv_adj_back_faces_ids;
    std::map<int,int > recv_adj_back_Nrhos;
    std::map<int,double* > recv_adj_back_rhos;
    int n_adj_vert_recv_back;
    int n_adj_face_recv_back;
    
    // This sends the right vertices of the requested elements to correct processor.
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = send_adj_verts_IDs.begin(); it != send_adj_verts_IDs.end(); it++)
            {
                int nv_adj_send       = it->second.size();
                int dest = it->first;
                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 98760000+1000*dest, comm);
                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 19999*9876+dest*8888,comm);
            
                int nf_adj_send = send_adj_faces_IDs[it->first].size();
                MPI_Send(&nf_adj_send, 1, MPI_INT, dest, 3333*9876+dest*8888,comm);
                MPI_Send(&send_adj_faces_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 2222*9876+dest*8888,comm);
                int n_adj_rhos = send_adj_rhos[it->first].size();
                MPI_Send(&n_adj_rhos, 1, MPI_INT, dest, 4444*9876+dest*8888,comm);
                MPI_Send(&send_adj_rhos[it->first][0], n_adj_rhos, MPI_DOUBLE, dest, 5555*9876+dest*8888,comm);

            }
        }
        if(sobj_el->RecvRankFromRank[q].find( rank ) != sobj_el->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 98760000+1000*rank, comm, MPI_STATUS_IGNORE);
            int* recv_adj_back_arr_ids = new int[n_adj_vert_recv_back];
            MPI_Recv(&recv_adj_back_arr_ids[0], n_adj_vert_recv_back, MPI_INT, q, 19999*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            int n_adj_face_recv_back;
            MPI_Recv(&n_adj_face_recv_back, 1, MPI_INT, q, 3333*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            int* recv_adj_back_arr_face_ids = new int[n_adj_face_recv_back];
            MPI_Recv(&recv_adj_back_arr_face_ids[0], n_adj_face_recv_back, MPI_INT, q, 2222*9876+rank*8888, comm,   MPI_STATUS_IGNORE);
            
            int n_adj_rho_recv_back;
            MPI_Recv(&n_adj_rho_recv_back, 1, MPI_INT, q, 4444*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            double* recv_adj_back_arr_rho = new double[n_adj_rho_recv_back];
            MPI_Recv(&recv_adj_back_arr_rho[0], n_adj_rho_recv_back, MPI_DOUBLE, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);

             
            recv_adj_back_Nverts[q]     = n_adj_vert_recv_back;
            recv_adj_back_verts_ids[q]  = recv_adj_back_arr_ids;
            
            recv_adj_back_Nfaces[q]     = n_adj_face_recv_back;
            recv_adj_back_faces_ids[q]  = recv_adj_back_arr_face_ids;
            
            recv_adj_back_Nrhos[q]      = n_adj_rho_recv_back;
            recv_adj_back_rhos[q]       = recv_adj_back_arr_rho;
    
        }
    }
        
    int TotNvert_adj_recv = 0;
    int TotNface_adj_recv = 0;
    int TotNrho_adj_recv  = 0;
    std::map<int,int >::iterator itm;
    std::vector<int> adj_verts;
    for(itm=recv_adj_back_Nverts.begin();itm!=recv_adj_back_Nverts.end();itm++)
    {
        TotNvert_adj_recv = TotNvert_adj_recv+itm->second;
        for(int i=0;i<itm->second;i++)
        {
            adj_verts.push_back(recv_adj_back_verts_ids[itm->first][i]);
        }
    }
    
    std::vector<int> adj_faces;
    for(itm=recv_adj_back_Nfaces.begin();itm!=recv_adj_back_Nfaces.end();itm++)
    {
        TotNface_adj_recv = TotNface_adj_recv+itm->second;
        for(int i=0;i<itm->second;i++)
        {
            adj_faces.push_back(recv_adj_back_faces_ids[itm->first][i]);
        }
    }
    
    std::vector<double> adj_rhos;
    for(itm=recv_adj_back_Nrhos.begin();itm!=recv_adj_back_Nrhos.end();itm++)
    {
        TotNrho_adj_recv = TotNrho_adj_recv+itm->second;
        for(int i=0;i<itm->second;i++)
        {
            adj_rhos.push_back(recv_adj_back_rhos[itm->first][i]);
        }
    }
    
    //std::cout << "TotNelem_adj_recv " << TotNelem_adj_recv << " should be equal to " << TotNrho_adj_recv << std::endl;
    int Nel_extra2 = itel; 
    //std::cout << " Compare " <<  TotNelem_adj_recv << " " << itel << " " << adj_rhos.size() << std::endl;
    int cnt_v_adj = 0;
    int cnt_f_adj = 0;
    for(int i=0;i<itel;i++)
    {
        
        for(int k=0;k<8;k++)
        {
            int v_id_n = adj_verts[cnt_v_adj];
            r = FindRank(new_offsets,size,v_id_n);
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id_n);
                //part_v.push_back(r);
                
                if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
                {
                    rank2req_vert[r].push_back(v_id_n); // add vertex to rank2req_vert map.
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n); // add the vertex to list that is already available on rank.
                    vloc++;
                }
            }
            cnt_v_adj++;
            
            if(k<6)
            {
                int f_id_n = adj_faces[cnt_f_adj];
                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
                {
                    unique_faceIDs_on_rank_set.insert(f_id);
                    faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                    floc++;
                }
                cnt_f_adj++;
            }
        }
        //part_elem2verts.push_back(elem);
        //elem_set.insert(TotAdj_El_IDs[i]);
    }
    
    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    
    
    
    // ==========================================================================================
    // ==========================================================================================
    // ==========================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // ==========================================================================================
    // ==========================================================================================
    // ==========================================================================================
    
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* sobj = DoScheduling(rank2req_vert,comm);
    
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_vert.begin(); it != rank2req_vert.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                
                int destination = dest;
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                
                i++;
            }
        }
        else if (sobj->SendFromRank2Rank[q].find( rank ) != sobj->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nverts;
    std::map<int,double* > recv_back_verts;
    std::map<int,int* > recv_back_verts_ids;
    int n_recv_back;
    //double* recv_back_arr = new double[10];
    
    
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                double* vert_send = new double[nv_send*3];
                offset_xcn        = xcn_parstate->getOffset(rank);
                nloc_xcn          = xcn_parstate->getNloc(rank);
                for(int u=0;u<it->second.size();u++)
                {
		    if(it->second[u]>(offset_xcn+nloc_xcn))
	       	    {
			std::cout << "Out of bounds > " << xcn->getNglob() << " " << it->second[u] << " " << offset_xcn+nloc_xcn << " " << offset_xcn << " " << nloc_xcn << std::endl;
		    }
		    else if(it->second[u]< offset_xcn)
		    {

			std::cout << "Out of bounds < " << xcn->getNglob() << " " << it->second[u] << " "      << offset_xcn+nloc_xcn << " " << offset_xcn << " " << nloc_xcn << std::endl;
		    } 
                    vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                    vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                    vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);
                }
                
                int dest = it->first;
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&vert_send[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, comm);
                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
		delete[] vert_send;
            }
        }
        if(sobj->RecvRankFromRank[q].find( rank ) != sobj->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, comm, MPI_STATUS_IGNORE);
            
            double* recv_back_arr = new double[n_recv_back*3];
            int* recv_back_arr_ids = new int[n_recv_back];
            //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr_ids[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q]     = n_recv_back;
            recv_back_verts[q]      = recv_back_arr;
            recv_back_verts_ids[q]  = recv_back_arr_ids;
        
	}
    }
   
    int vfor = 0;
    std::map<int,double* >::iterator it_f;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int c  = 0;
        vfor=vfor+recv_back_Nverts[it_f->first];

    }
    double* part_verts_arr = new double[3*(vloc+vfor)];

    int gvid=0;
    int lvid=0;
    
    for(int m=0;m<vloc;m++)
    {
        gvid = vertIDs_on_rank[m];
       
        part_verts_arr[m*3+0] = xcn->getVal(gvid-xcn_o,0);
        part_verts_arr[m*3+1] = xcn->getVal(gvid-xcn_o,1);
        part_verts_arr[m*3+2] = xcn->getVal(gvid-xcn_o,2);
        Vert V;
        
        V.x = xcn->getVal(gvid-xcn_o,0);
        V.y = xcn->getVal(gvid-xcn_o,1);
        V.z = xcn->getVal(gvid-xcn_o,2);
        
        LocalVerts.push_back(V);
        LocalVert2GlobalVert[lvid] = gvid;
        GlobalVert2LocalVert[gvid] = lvid;
        lvid++;
    }
    
    int o = 3*vloc;
    int m = 0;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = rank2req_vert[it_f->first][u];
            
            
            part_verts_arr[o+m*3+0] = it_f->second[u*3+0];
            part_verts_arr[o+m*3+1] = it_f->second[u*3+1];
            part_verts_arr[o+m*3+2] = it_f->second[u*3+2];
            Vert V;
            
            V.x = it_f->second[u*3+0];
            V.y = it_f->second[u*3+1];
            V.z = it_f->second[u*3+2];
            
            LocalVerts.push_back(V);
            
            LocalVert2GlobalVert[lvid]=gvid;
            GlobalVert2LocalVert[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }
    
    NlocVerts = LocalVerts.size();
    // ================================== Faces on Rank =========================================
    
    int lfid = 0;
    int gfid = 0;
    for(int m=0;m<floc;m++)
    {
        gfid = faceIDs_on_rank[m];
    
        LocalFace2GlobalFace[lfid] = gfid;
        GlobalFace2LocalFace[gfid] = lfid;
        lfid++;
    }
    
    // ================================== Faces on Rank =========================================
    NlocElem             = loc_elem.size()+Nel_extra+Nel_extra2;
    LocalElem2GlobalVert = new Array<int>(NlocElem,8);
    LocalElem2LocalVert  = new Array<int>(NlocElem,8);
    
    U0Elem               = new Array<double>(NlocElem,1);
    U0Vert               = new Array<double>(LocalVerts.size(),1);
    ElemPart             = new Array<int>(NlocElem,1);

    int glob_v = 0;
    int loc_v  = 0;
    int glob_f = 0;
    int loc_f  = 0;
    double rho_v = 0;
    int loc_el = 0;
    std::map<int,std::vector<double> > collect_var;
    for(int m=0;m<loc_elem.size();m++)
    {
        el_id = loc_elem[m];
        rho_v = loc_rho[m];
        U0Elem->setVal(m,0,rho_v);
        ElemPart->setVal(m,0,el_id);
        LocalElement2GlobalElement[loc_el] = el_id;
        GlobalElement2LocalElement[el_id] = loc_el;
        loc_el++;
        
        for(int p=0;p<8;p++)
        {
            //GlobalFace2LocalFace
            //LocalVert2GlobalVert
            glob_v = ien->getVal(el_id-ien_o,p);
            loc_v = GlobalVert2LocalVert[glob_v];
            LocalElem2GlobalVert->setVal(m,p,glob_v);
            LocalElem2LocalVert->setVal(m,p,loc_v);
            collect_var[loc_v].push_back(rho_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            
//            if(rank == 0)
//            {
//            std::cout << glob_v << " ";
//            }
            if(p<6)
            {
                glob_f = ief->getVal(el_id-ien_o,p);
                loc_f  = GlobalFace2LocalFace[glob_f];
                
                Elem2GlobalFace[el_id].push_back(glob_f);
                GlobalFace2Elem[glob_f].push_back(el_id);
            }
        }
//        if(rank == 0)
//        {
//        std::cout << std::endl;
//        }

    }
    o = loc_elem.size();
    int cnv = 0;
    int cnf = 0;
    for(int m=0;m<Nel_extra;m++)
    {
        el_id = TotRecvElement_IDs[m];
        LocalElement2GlobalElement[loc_el] = el_id;
        GlobalElement2LocalElement[el_id] = loc_el;
        loc_el++;
        rho_v = TotRecvElement_rhos[m];
        U0Elem->setVal(m+o,0,rho_v);
        ElemPart->setVal(m+o,0,el_id);

        for(int p=0;p<8;p++)
        {
            glob_v = TotRecvElement_IDs_v[cnv];
            loc_v = GlobalVert2LocalVert[glob_v];
            LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            LocalElem2LocalVert->setVal(m+o,p,loc_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            collect_var[loc_v].push_back(rho_v);
            cnv++;
            if(p<6)
            {
                glob_f = TotRecvElement_IDs_f[cnf];
                loc_f  = GlobalFace2LocalFace[glob_f];

                Elem2GlobalFace[el_id].push_back(glob_f);
                GlobalFace2Elem[glob_f].push_back(el_id);
                cnf++;
            }
        }
    }
    
    o = loc_elem.size()+Nel_extra;
    cnv = 0;
    cnf = 0;
    int ofs = xcn_parstate->getOffset(rank);
    int nlo = xcn_parstate->getNloc(rank);

    for(int m=0;m<Nel_extra2;m++)
    {
        el_id = adj_elements[m];
        LocalElement2GlobalElement[loc_el] = el_id;
        GlobalElement2LocalElement[el_id] = loc_el;
        loc_el++;
        rho_v = adj_rhos[m];
        U0Elem->setVal(m+o,0,rho_v);
        ElemPart->setVal(m+o,0,el_id);

        for(int p=0;p<8;p++)
        {
            glob_v = adj_verts[cnv];
            loc_v  = GlobalVert2LocalVert[glob_v];
            
            LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            LocalElem2LocalVert->setVal(m+o,p,loc_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            collect_var[loc_v].push_back(rho_v);
            cnv++;
            if(p<6)
            {
                glob_f = adj_faces[cnf];
                loc_f  = GlobalFace2LocalFace[glob_f];

                Elem2GlobalFace[el_id].push_back(glob_f);
                GlobalFace2Elem[glob_f].push_back(el_id);
                cnf++;
            }
        }
    }
    
    std::map<int,std::vector<double> >::iterator it_rhos;
    double sum = 0;
    int c = 0;
    
    for(it_rhos=collect_var.begin();it_rhos!=collect_var.end();it_rhos++)
    {
        sum = 0;
        for(int q = 0;q<it_rhos->second.size();q++)
        {
            sum = sum + it_rhos->second[q];
        }
        U0Vert->setVal(c,0,sum/it_rhos->second.size());
        c++;
    }
    
    //std::cout << rank << " Partitioning on Rank = " << (double) not_on_rank/on_rank << std::endl;
}

inline int Partition::getNlocElem()
{
    return NlocElem;
}

inline int Partition::getNlocVerts()
{
    return NlocVerts;
}

inline int* Partition::getXadj()
{
    return xadj;
}
inline int* Partition::getAdjcny()
{
    return adjcny;
}
inline ParArray<int>* Partition::getPart()
{
    return part;
}
inline std::vector<Vert> Partition::getLocalVerts()
{
    return LocalVerts;
}
inline Vert Partition::getLocalVert(int v_loc_id)
{
    return LocalVerts[v_loc_id];
}
inline Array<int>* Partition::getLocalElem2GlobalVert()
{
    return LocalElem2GlobalVert;
}
inline Array<int>* Partition::getLocalElem2LocalVert()
{
    return LocalElem2LocalVert;
}
inline std::map<int,int> Partition::getLocalVert2GlobalVert()
{
    return LocalVert2GlobalVert;
}
inline std::map<int,int> Partition::getGlobalVert2LocalVert()
{
    return GlobalVert2LocalVert;
}


inline Array<int>* Partition::getLocalElem2GlobalFace()
{
    return LocalElem2GlobalFace;
}
inline Array<int>* Partition::getLocalElem2LocalFace()
{
    return LocalElem2LocalFace;
}
inline std::map<int,int> Partition::getLocalFace2GlobalFace()
{
    return LocalFace2GlobalFace;
}
inline std::map<int,int> Partition::getGlobalFace2LocalFace()
{
    return GlobalFace2LocalFace;
}
inline std::map<int, std::vector<int> > Partition::getElem2GlobalFace()
{
    return Elem2GlobalFace;
}
inline std::map<int, std::vector<int> > Partition::getGlobalFace2Elem()
{
    return GlobalFace2Elem;
}
inline std::set<int> Partition::getElemSet()
{
    return elem_set;
}
inline Array<double>* Partition::getUelem()
{
    return U0Elem;
}
inline double Partition::getU0atGlobalElem(int gelem)
{
    int elem = GlobalElement2LocalElement[gelem];
    return U0Elem->getVal(elem,0);
}
inline Array<double>* Partition::getUvert()
{
    return U0Vert;
}
inline std::map<int,std::vector<int> > Partition::getGlobElem2LocVerts()
{
    return globElem2locVerts;
}
inline std::map<int,std::vector<int> > Partition::getGlobElem2GlobVerts()
{
    return globElem2globVerts;
}
inline std::map<int,int> Partition::getGlobalElement2LocalElement()
{
    return GlobalElement2LocalElement;
}
inline std::map<int,int> Partition::getLocalElement2GlobalElement()
{
    return LocalElement2GlobalElement;
}

//ParallelState* getXcnParallelState();
//ParallelState* getIenParallelState();
//ParallelState_Parmetis* getParallelStateParmetis();
#endif
