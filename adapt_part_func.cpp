#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#include <iostream>
#include "adapt_part_func.h"

using namespace std;






GathervObject* GetGathervObject(int nloc, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int* locs     = new int[size];
    int* reduced_locs = new int[size];

    for(int i=0;i<size;i++)
    {
        reduced_locs[i]  = 0;
        
        if(i==rank)
        {
            locs[i]  = nloc;
        }
        else
        {
            locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(locs, reduced_locs, size, MPI_INT, MPI_SUM, comm);
    
    int* reduced_offsets = new int[size];
    reduced_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        reduced_offsets[i+1]=reduced_offsets[i]+reduced_locs[i];
    }
    
    int length = reduced_offsets[size-1]+reduced_locs[size-1];
    
    GathervObject* gObj = new GathervObject;
    gObj->data = new int[length];
    gObj->nlocs = reduced_locs;
    gObj->offsets = reduced_offsets;
    gObj->size    = size;
    gObj->length  = length;
    
    return gObj;
    
}


std::vector<int> GetAdjacencyForUS3D_V4(ParArray<int>* ief, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = ief->getNrow();
    int ncol = ief->getNcol();

//    Array<int>*type;
    std::vector<int> ief_copy(nrow*(ncol-1));
    set<int> unique_faces;
    set<int> double_faces;
    std::vector<int> d_faces;
    std::vector<int> u_faces;
    
   
    int fid;
    
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            fid = fabs(ief->getVal(i,j+1))-1;
            
            if(unique_faces.find(fid) == unique_faces.end())
            {
                unique_faces.insert(fid);
                u_faces.push_back(fid);
            }
            else
            {
                double_faces.insert(fid);
                d_faces.push_back(fid);
            }
        }
    }

    
    int exter_size = u_faces.size();

    int* exter_nlocs                 = new int[size];
    int* exter_offset                = new int[size];
    int* reduced_exter_nlocs             = new int[size];
    int* reduced_exter_offset            = new int[size];
    

    
    for(int i=0;i<size;i++)
    {
        reduced_exter_nlocs[i]  = 0;
        reduced_exter_offset[i] = 0;
        if(i==rank)
        {
            exter_nlocs[i] = exter_size;
        }
        else
        {
            exter_nlocs[i]  = 0;
            exter_offset[i] = 0;
        }
    }
    
    MPI_Allreduce(exter_nlocs,
                  reduced_exter_nlocs,
                  size,
                  MPI_INT,
                  MPI_SUM,
                  comm);
    
    
    reduced_exter_offset[0]=0;
        
    
    //
    for(int i=0;i<size-1;i++)
    {
        reduced_exter_offset[i+1]=reduced_exter_offset[i]+reduced_exter_nlocs[i];
    }
    //
    
    int nexter_tot = reduced_exter_offset[size-1]+reduced_exter_nlocs[size-1];
    
    /*
    std::vector<int> recv(nexter_tot);
     
    MPI_Allgatherv(&u_faces[0],
                   exter_size,
                   MPI_INT,
                   &recv[0],
                   reduced_exter_nlocs,
                   reduced_exter_offset,
                   MPI_INT, comm);
    */
    

    //start = std::clock();
    
    int* uf_arr = new int[u_faces.size()];
    for(int i=0;i<u_faces.size();i++)
    {
        uf_arr[i] = u_faces[i];
    }
    
    std::vector<int> recv2 = FindDuplicatesInParallel(uf_arr, u_faces.size(), nexter_tot, comm);
    
    
    
    delete[] uf_arr;
    delete[] exter_nlocs;
    delete[] exter_offset;
    delete[] reduced_exter_nlocs;
    delete[] reduced_exter_offset;
    
    
    return recv2;


}




ParArrayOnRoot* GatherVecToRoot(std::vector<int> locvec, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int* locs     = new int[size];
    int* reduced_locs = new int[size];

    for(int i=0;i<size;i++)
    {
        reduced_locs[i]  = 0;
        
        if(i==rank)
        {
            locs[i]  = locvec.size();
        }
        else
        {
            locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(locs, reduced_locs, size, MPI_INT, MPI_SUM, comm);
    
    int* reduced_offsets = new int[size];
    reduced_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        reduced_offsets[i+1]=reduced_offsets[i]+reduced_locs[i];
    }
    
    int tot = reduced_offsets[size-1]+reduced_locs[size-1];
    
    ParArrayOnRoot* parr_root = new ParArrayOnRoot;
    parr_root->data           = new int[tot];
    parr_root->nlocs          = reduced_locs;
    parr_root->offsets        = reduced_offsets;
    parr_root->size           = size;
    parr_root->length         = tot;
    
    MPI_Gatherv(&locvec[0],
                locvec.size(),
                MPI_INT,
                &parr_root->data[0],
                parr_root->nlocs,
                parr_root->offsets,
                MPI_INT,0, comm);
    
    return parr_root;
}
/*
int* GetPartitionInfo(ParArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm)
{
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    int nloc = nrow;
    
    int N = ien->getNglob();
    ParArray<int>* ien_copy = new ParArray<int>(N, ncol-1, comm);

    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            ien_copy->setVal(i,j,ien->getVal(i,j+1)-1);
        }
    }

    //=================================================================
    //=================================================================
    //=================================================================
    
    ParVar_ParMetis* pv_parmetis = CreateParallelDataParmetis(ien_copy,comm,8);
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    

    
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt  = NULL;
    int np           = size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    int* part = new int[nloc];

    ParMETIS_V3_PartMeshKway(pv_parmetis->elmdist,
                             pv_parmetis->eptr,
                             pv_parmetis->eind,
                             elmwgt, wgtflag, numflag,
                             ncon, ncommonnodes, nparts,
                             tpwgts, ubvec, options,
                             &edgecut, part, &comm);
    
    
    int L = pv_parmetis->elmdist[size];
    int* part_collect_on_root = NULL;
    
    if(rank == 0)
    {
        part_collect_on_root=new int[L];
    }
    
    MPI_Gatherv(&part[0],
                nloc,
                MPI_INT,
                &part_collect_on_root[0],
                ien->getParallelState()->getNlocs(),
                ien->getParallelState()->getOffsets(),
                MPI_INT,0, comm);

    
    delete ien_copy;
    delete pv_parmetis;
    delete[] part;
    return part_collect_on_root;
}
*/
Partition_old* CollectVerticesPerRank(ParArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    
    //=================================================================
    //=================================================================
    //=================================================================
    ParallelState_Parmetis* pstate_parmetis = new ParallelState_Parmetis(ien,comm,8);
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;

    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj,&adjncy,&comm);
    
    //=================================================================
    //=================================================================
    //=================================================================
    
    std::map<int,int> loc2glob;
    std::map<int,int> glob2loc;
    Array<int>* loc2glob_vert = new Array<int>(nrow,ncol-1);
    Array<int>* glob2loc_vert = new Array<int>(nrow,ncol-1);

    std::vector<int> loc_elems;
    set<int> loc_elems_set;
    int gid = 0;
    int lid = 0;
    
    for(int i=0;i<nrow;i++)
    {
        int glob_id = i+pstate_parmetis->getElmdistAtRank(rank);
        if ( loc_elems_set.find( glob_id ) == loc_elems_set.end() )
        {
            loc_elems.push_back(glob_id);
            loc_elems_set.insert(glob_id);
            loc2glob[lid] = gid;
            glob2loc[gid] = lid;
            lid++;
        }
        for(int j=xadj[i];j<xadj[i+1];j++)
        {
            if ( loc_elems_set.find( adjncy[j] ) == loc_elems_set.end() )
            {
                loc_elems.push_back(adjncy[j]);
                loc_elems_set.insert(adjncy[j]);
                loc2glob[lid] = gid;
                glob2loc[gid] = lid;
                lid++;
            }
        }
    }
    
    // This vector is empty on all other procs except root;
    ParArrayOnRoot* gathereduced_on_root = GatherVecToRoot(loc_elems, comm);
    
    
    int* nlocs = new int[size];
    int* offset = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs[i] = gathereduced_on_root->nlocs[i]*3;
        offset[i] = gathereduced_on_root->offsets[i]*3;
    }
    
    Partition_old* parti = new Partition_old;
    
    parti->Verts = new Array<double>(loc_elems.size(),3);
    
    parti->loc2glob_Vmap = loc2glob;
    parti->glob2loc_Vmap = glob2loc;
    
    parti->loc2glob_Varr = loc2glob_vert;
    parti->glob2loc_Varr = glob2loc_vert;
    
    parti->xadj   = xadj;
    parti->adjncy = adjncy;
    parti->ien    = ien;
    parti->ndim   = 3;
    
    double * verts = NULL;
    if(rank == 0)
    {
        verts = new double[gathereduced_on_root->length*3];

        for(int i=0;i<gathereduced_on_root->length;i++)
        {
            verts[i*3+0] = xcn_r->getVal(gathereduced_on_root->data[i],0);
            verts[i*3+1] = xcn_r->getVal(gathereduced_on_root->data[i],1);
            verts[i*3+2] = xcn_r->getVal(gathereduced_on_root->data[i],2);
        }
    }
 
    MPI_Scatterv(&verts[0], nlocs, offset,
                 MPI_DOUBLE, &parti->Verts->data[0],
                 loc_elems.size()*3, MPI_DOUBLE, 0, comm);
    
    delete[] xadj;
    delete[] adjncy;
    delete pstate_parmetis;
    return parti;
}


ParArray<int>* DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, MPI_Comm comm)
{
    ParArray<int>*  part_arr;
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
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
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
    int* part = new int[nloc];
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL; 
    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj,&adjncy,&comm);
    
    /*
    for(int u=0;u<nloc;u++)
    {
    	part[u] = rank;
    }*/ 
    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(), 
	                       xadj, adjncy, 
			       elmwgt, adjwgt, 
		               vsize, wgtflag, 
			       numflag, ncon, nparts,
			       tpwgts, ubvec, itr, options, 
			       &edgecut, part, &comm);
 
    /*ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj,
                         adjncy,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part, &comm);
     */
    part_arr = new ParArray<int>(ien->getNglob(),1,comm);
    part_arr->data = part;
     
    //delete pstate_parmetis;
    
    return part_arr;
    
}

//================================================================================
// This function determines a map that gets the unique list of elements for that need to be requested from a given rank other than current rank.
PartitionStruct* DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* part, ParArray<double>* xcn, ParArray<double>* variables, ParallelState* xcn_parstate, MPI_Comm comm)
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
    std::map<int,std::vector<int> > verts_to_send_to_ranks;   
    std::map<int,std::vector<double> > rho_to_send_to_ranks;

    set<int> unique_verts_on_rank_set;
    std::vector<int> unique_verts_on_rank_vec;
    
    set<int> u_verts_part_set;
    std::vector<int> u_verts_part_vec;
    std::map<int,std::vector<int> > rank2req_vert;

    std::vector<int> vert_on_rank;
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
    std::map<int,int> v_loc2glob;
    std::map<int,int> v_glob2loc;
    int vloc = 0;
    std::vector<int> loc_elem;
    std::vector<double> loc_rho;
    int ien_o = part->getOffset(rank);
    double rho = 0.0;
    for(i=0;i<part->getNrow();i++)
    {
        p_id  = part->getVal(i,0);
        el_id = part->getOffset(rank)+i;
        rho = variables->getVal(i,0);
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            rho_to_send_to_ranks[p_id].push_back(rho);
            
            for(int k=0;k<8;k++)
            {
                v_id = ien->getVal(i,k);
                verts_to_send_to_ranks[p_id].push_back(v_id);
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;
            for(int k=0;k<8;k++)// looping over the vertices for element "i".
            {
                v_id = ien->getVal(i,k);
                //elem.push_back(v_id);
                if(unique_verts_on_rank_set.find( v_id ) == unique_verts_on_rank_set.end())// find the unique vertices that need to be send to other partitions.
                {
                    unique_verts_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(xcn_parstate->getOffsets(),size,v_id);

                    if (r!=rank)// if vertex is present on other rank, add it to vert_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        vert_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vloc++;
                    }
                    lv_id++;
                }
            }
            loc_elem.push_back(el_id);
            loc_rho.push_back(rho);
        }
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
        sendNelemFromRank2Rank_Global[i] = 0;
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
    
    
    std::map<int, set<int> > sendFromRank2Rank_e;
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
    
    int* TotRecvElement_IDs   = new int[TotNelem_recv];
    double* TotRecvElement_rhos   = new double[TotNelem_recv];

    int* TotRecvElement_IDs_v = new int[TotNelem_recv*8];
    
    for(int i =0;i<TotNelem_recv;i++)
    {
        TotRecvElement_IDs[i] = 0;
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
                int dest            = it->first;
                                
                MPI_Send(&n_req, 1, MPI_INT, dest, dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                MPI_Send(&rho_to_send_to_ranks[it->first][0], n_req, MPI_DOUBLE, dest, 20000+100+dest*2, comm);
                
                MPI_Send(&n_req_v, 1, MPI_INT, dest, 9000+dest, comm);
                MPI_Send(&verts_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                
                i++;
            }
        }
        else if (sendFromRank2Rank_e[q].find( rank ) != sendFromRank2Rank_e[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_IDs[RecvAlloc_offset_map_e[q]], n_req_recv, MPI_INT, q, 100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_rhos[RecvAlloc_offset_map_e[q]], n_req_recv, MPI_DOUBLE, q, 20000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, 9000+rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&TotRecvElement_IDs_v[RecvAlloc_offset_map_e[q]*8], n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
        }
    }

    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    int Nel_extra = TotNelem_recv;
    int cnt = 0;
    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;
        for(int k=0;k<8;k++)
        {
            int v_id_n = TotRecvElement_IDs_v[cnt];
            elem.push_back(v_id_n);
            r = FindRank(xcn_parstate->getOffsets(),size,v_id_n);
            
            if(unique_verts_on_rank_set.find( v_id_n ) == unique_verts_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_verts_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id_n);
                //part_v.push_back(r);
                
                if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vert_on_rank map.
                {
                    rank2req_vert[r].push_back(v_id_n); // add vertex to rank2req_vert map.
                }
                else
                {
                    vert_on_rank.push_back(v_id_n); // add the vertex to list that is already available on rank.
                    vloc++;
                }
            }
            cnt++;
        }
        //part_elem2verts.push_back(elem);
        //loc_elem.push_back(TotRecvElement_IDs[i]);
    }
    
    // ==========================================================================================================
    // ==========================================================================================================
    // ==========================================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // ==========================================================================================================
    // ==========================================================================================================
    // ==========================================================================================================
    
    
    int nRank_reqVerts          = rank2req_vert.size(); // The number of ranks from which current rank requests vertices.
    int* reduced_nRank_reqVerts = new int[size]; // Defined memory for a reduced array so that all ranks are going to be aware of which information is required from each other.
    int* arr_nRank_reqVerts     = new int[size]; // Defining memory to store local requesting information.
    
    for(i=0;i<size;i++)
    {
        reduced_nRank_reqVerts[i] = 0;
        
        if(i==rank)
        {
            arr_nRank_reqVerts[i] = nRank_reqVerts+1;
        }
        else
        {
            arr_nRank_reqVerts[i] = 0;
        }
    }
    
    
    // Reduce the requesting schedule for the vertices to all processors.
    MPI_Allreduce(arr_nRank_reqVerts,
                  reduced_nRank_reqVerts,
                  size, MPI_INT, MPI_SUM, comm);
    
    
    int* reduced_nRank_reqVerts_offset = new int[size];// Define an offset array for the schedule in order to be able to gather the local schedules for each rank to a global schedule array.
    
    offset = 0;
    for(i=0;i<size;i++)
    {
        reduced_nRank_reqVerts_offset[i] = offset;
        offset = offset+reduced_nRank_reqVerts[i];
    }
    int nTot_reqVerts = 0;
    int nRank_reqVerts_p_one = nRank_reqVerts+1; // This size is added by one since we add the rank number to the array.
    
    MPI_Allreduce(&nRank_reqVerts_p_one, &nTot_reqVerts, 1, MPI_INT, MPI_SUM, comm);// Determine the total length of the "schedule" array.
    int* sendFromRank2Rank_vert_Global = new int[nTot_reqVerts]; // This array is laid out as follows:
    // first the fromRank is noted which is followed by the IDs for the several ranks FromRank is sending to.
    int* sendNvertFromRank2Rank_Global = new int[nTot_reqVerts]; // This array is layout as follows:
    // first the fromRank is noted which is followed by the Nverts is listed for the several ranks FromRank is sending to.
    
    for(i=0;i<nTot_reqVerts;i++)
    {
        sendFromRank2Rank_vert_Global[i] = 0;
        sendNvertFromRank2Rank_Global[i] = 0;
    }
    
    int* reqRank_fromRank_Vert   = new int[nRank_reqVerts+1];
    int* reqNverts_from_rank     = new int[nRank_reqVerts+1];
    reqRank_fromRank_Vert[0]     = rank;
    reqNverts_from_rank[0]       = -1;
    t = 1;
    std::map<int,std::vector<int> >::iterator it2;
    for(it2=rank2req_vert.begin();it2!=rank2req_vert.end();it2++)
    {
        reqRank_fromRank_Vert[t]     = it2->first;
        reqNverts_from_rank[t]       = it2->second.size();
        t++;
    }

    MPI_Allgatherv(&reqRank_fromRank_Vert[0], nRank_reqVerts_p_one, MPI_INT,
                   &sendFromRank2Rank_vert_Global[0],
                   reduced_nRank_reqVerts, reduced_nRank_reqVerts_offset, MPI_INT,comm);
    
    MPI_Allgatherv(&reqNverts_from_rank[0], nRank_reqVerts_p_one, MPI_INT,
                   &sendNvertFromRank2Rank_Global[0],
                   reduced_nRank_reqVerts, reduced_nRank_reqVerts_offset, MPI_INT,comm);
    

    std::map<int, set<int> > sendFromRank2Rank_v_set;
    std::map<int, std::vector<int> > recvRankFromRank_map_v_vec;
    std::map<int, std::vector<int> > recvRankFromRank_map_Nvert_v_vec;
    std::map<int, set<int> > recvRankFromRank_map_v_set;
    //=========================================
    for(i=0;i<size;i++)
    {
        int of = reduced_nRank_reqVerts_offset[i];
        int nl = reduced_nRank_reqVerts[i];
        
        for(int j=of+1;j<of+nl;j++)
        {
            sendFromRank2Rank_v_set[sendFromRank2Rank_vert_Global[of]].insert(sendFromRank2Rank_vert_Global[j]);
            
            recvRankFromRank_map_v_set[sendFromRank2Rank_vert_Global[j]].insert(sendFromRank2Rank_vert_Global[of]);
            recvRankFromRank_map_v_vec[sendFromRank2Rank_vert_Global[j]].push_back(sendFromRank2Rank_vert_Global[of]);
            recvRankFromRank_map_Nvert_v_vec[sendFromRank2Rank_vert_Global[j]].push_back(sendNvertFromRank2Rank_Global[j]);
        }
    }
    
    std::map<int,std::map<int,int> > glob_rank_rank2req_vert; // This map links current rank to rank2vert map.
    
    for(it=recvRankFromRank_map_v_vec.begin();it!=recvRankFromRank_map_v_vec.end();it++)
    {
        //it->first is the value of the receiving rank
        //it->second are the ids of the ranks that send data to the receiving rank.
        for(int k=0;k<it->second.size();k++)
        {
            glob_rank_rank2req_vert[it->first][recvRankFromRank_map_v_vec[it->first][k]] = recvRankFromRank_map_Nvert_v_vec[it->first][k];
        }
    }
    
    std::map<int,int> local_rank_rank2req_vert = glob_rank_rank2req_vert[rank];
    
    std::map<int, std::vector<int> >::iterator it8;

    
    int* recv_offset2 = new int[local_rank_rank2req_vert.size()+1];
    recv_offset2[0]   = 0;
    int* recv_loc2    = new int[local_rank_rank2req_vert.size()];
    int TotNvert_recv    = 0;
    
    std::map< int, int> RecvAlloc_offset_map_v;
    std::map< int, int> loc_map2;
    std::map< int, int>::iterator it_loc2;
    i = 0;
    int offs = 0;
    
    for(it_loc2=local_rank_rank2req_vert.begin();it_loc2!=local_rank_rank2req_vert.end();it_loc2++)
    {
        recv_loc2[i]                = it_loc2->second;
        recv_offset2[i+1]           = recv_offset2[i]+recv_loc2[i];
        loc_map2[it_loc2->first]    = recv_loc2[i];
        RecvAlloc_offset_map_v[it_loc2->first] = recv_offset2[i];
        offs                        = offs+it_loc2->second;
        TotNvert_recv               = TotNvert_recv+it_loc2->second;
        i++;
    }
    
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    int* TotRecvVert_IDs = new int[TotNvert_recv];
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
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
                
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                
                i++;
            }
        }
        else if (sendFromRank2Rank_v_set[q].find( rank ) != sendFromRank2Rank_v_set[q].end())
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
    std::map<int,int  > recv_back_Nverts;
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
                offset_xcn        = xcn->getOffset(rank);
                nloc_xcn          = xcn->getNloc(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                    vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                    vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);
                }
                
                int dest = it->first;
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&vert_send[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, comm);
                //MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                delete[] vert_send;
            }
        }
        if(recvRankFromRank_map_v_set[q].find( rank ) != recvRankFromRank_map_v_set[q].end())
        {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, comm, MPI_STATUS_IGNORE);
            
            double* recv_back_arr = new double[n_recv_back*3];
            int* recv_back_arr_ids = new int[n_recv_back];
            //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&recv_back_arr_ids[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q] = n_recv_back;
            recv_back_verts[q]  = recv_back_arr;
            recv_back_verts_ids[q]  = recv_back_arr_ids;
        }
    }
    
    int vfor = 0;
    std::map<int,double* >::iterator it_f;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int c  = 0;
        vfor=vfor+recv_back_Nverts[it_f->first];
        //for(int u=0;u<Nv;u++)
        //{
            //v_id = rank2req_vert[it_f->first][u];
            //vert_on_rank.push_back(v_id);
            
            //V.x = it_f->second[u*3+0];
            //V.y = it_f->second[u*3+1];
            //V.z = it_f->second[u*3+2];
            
            //part_verts.push_back(V);
            //v_loc2glob[vloc]=v_id;
            //v_glob2loc[v_id]=vloc;
 	    //vfor++;
            //vloc++;
        //}
    }
    double* part_verts_arr = new double[3*(vloc+vfor)];
    //int* v_loc2glob = new int[vloc+vfor];
    //int* v_glob2loc = new int[vloc+vfor];
    int vid=0;
    int lid=0;
    std::vector<Vert> Verts;
    for(int m=0;m<vloc;m++)
    {  
        vid = vert_on_rank[m];
       
        part_verts_arr[m*3+0] = xcn->getVal(vid-xcn_o,0);
        part_verts_arr[m*3+1] = xcn->getVal(vid-xcn_o,1);
        part_verts_arr[m*3+2] = xcn->getVal(vid-xcn_o,2);
        Vert V;
        
        V.x = xcn->getVal(vid-xcn_o,0);
        V.y = xcn->getVal(vid-xcn_o,1);
        V.z = xcn->getVal(vid-xcn_o,2);
        
        Verts.push_back(V);
        v_loc2glob[lid] = vid;
        v_glob2loc[vid] = lid;
        lid++;
    }
    

    int o = 3*vloc;
    int m = 0;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {   
        int Nv = recv_back_Nverts[it_f->first];
        for(int u=0;u<Nv;u++)
        {  
            vid = rank2req_vert[it_f->first][u];
     
            part_verts_arr[o+m*3+0] = it_f->second[u*3+0];
            part_verts_arr[o+m*3+1] = it_f->second[u*3+1];
            part_verts_arr[o+m*3+2] = it_f->second[u*3+2];
            Vert V;
            
            V.x = it_f->second[u*3+0];
            V.y = it_f->second[u*3+1];
            V.z = it_f->second[u*3+2];
            
            Verts.push_back(V);
            
            v_loc2glob[lid]=vid;
            v_glob2loc[vid]=lid;
            
            m++;
            lid++;
            //vfor++;
            //vloc++;
        }
        
    }
    
    int Nloc_elem = loc_elem.size()+Nel_extra;
    Array<int>* part_El2Vert_glob = new Array<int>(Nloc_elem,8);
    Array<int>* part_El2Vert_loc  = new Array<int>(Nloc_elem,8);
    Array<double>* rho_part       = new Array<double>(Nloc_elem,1);
    Array<double>* rho_part_v     = new Array<double>(Nloc_elem,1);
    Array<int>* elem_part         = new Array<int>(Nloc_elem,1);

    int glob = 0;
    int loc = 0;
    double rho_v = 0;
    std::map<int,std::vector<double> > collect_var;
    for(int m=0;m<loc_elem.size();m++)
    {
        el_id = loc_elem[m];
        rho_v = loc_rho[m];
        rho_part->setVal(m,0,rho_v);
        elem_part->setVal(m,0,el_id);
        for(int p=0;p<8;p++)
        {
            glob = ien->getVal(el_id-ien_o,p);
            loc = v_glob2loc[glob];
            part_El2Vert_glob->setVal(m,p,glob);
            part_El2Vert_loc->setVal(m,p,loc);
            collect_var[loc].push_back(rho_v);
            
        }
    }
    o = loc_elem.size();
    int cn = 0;
    for(int m=0;m<Nel_extra;m++)
    {	
        el_id = TotRecvElement_IDs[m];
        rho_v = TotRecvElement_IDs[m];
        rho_part->setVal(m+o,0,rho_v);
        elem_part->setVal(m+o,0,el_id);
        
        for(int p=0;p<8;p++)
        {    
            glob = TotRecvElement_IDs_v[cn];
            loc = v_glob2loc[glob];
            part_El2Vert_glob->setVal(m+o,p,glob);
            part_El2Vert_loc->setVal(m+o,p,loc);
            collect_var[loc].push_back(rho_v);
            cn++;
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
        rho_part_v->setVal(c,0,sum/it_rhos->second.size());
        c++;
    }

    PartitionStruct* P               = new PartitionStruct;
    P->Verts                   = Verts;//part_verts;
    P->loc_elem2verts_glob     = part_El2Vert_glob;
    P->loc_elem2verts_loc      = part_El2Vert_loc;
    P->v_loc2glob              = v_loc2glob;
    P->v_glob2loc              = v_glob2loc;
    P->rho_elem                = rho_part; // This is the reduced average for each vert based on the related cell value.
    P->rho_vert                = rho_part_v; // This is the reduced average for each vert based on the related cell value.
    int tot_num_elem = 0;
    MPI_Allreduce(&Nloc_elem, &tot_num_elem, 1, MPI_INT, MPI_SUM, comm);
    
    if (rank == 0)
    {
        std::cout << "the total number of elements := " << tot_num_elem << " " << ien->getNglob() << std::endl;
    }
    
    
    
    /*

//
    
    
    // ==========================================================================================================
    // ==========================================================================================================
    // ==========================================================================================================
    
    
    
    
    
    // ---> First attempt by gathering the schedule on root and divide by root which requires long loop on root causing a serial bottleneck. <--- //
    
    ParArrayOnRoot* gathereduced_on_root = GatherVecToRoot(unique_verts_on_rank_vec, comm);
    
    int* nlocs_scat = new int[size];
    int* offset_scat = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs_scat[i]  = gathereduced_on_root->nlocs[i]*3;
        offset_scat[i] = gathereduced_on_root->offsets[i]*3;
    }
       
    double* verts = NULL;
    Array<double>* Verts_loc = new Array<double>(unique_verts_on_rank_vec.size(),3);
    if(rank == 0)
    {
        
        verts = new double[gathereduced_on_root->length*3];
        for(int i=0;i<gathereduced_on_root->length;i++)
        {
            verts[i*3+0] = xcn_on_root->getVal(gathereduced_on_root->data[i],0);
            verts[i*3+1] = xcn_on_root->getVal(gathereduced_on_root->data[i],1);
            verts[i*3+2] = xcn_on_root->getVal(gathereduced_on_root->data[i],2);
        }
        
    }
    Partition2* P2 = new Partition2;
    MPI_Scatterv(&verts[0], nlocs_scat, offset_scat, MPI_DOUBLE, &Verts_loc->data[0], unique_verts_on_rank_vec.size()*3, MPI_DOUBLE, 0, comm);
    */ 
    
    //delete[] verts;
    //delete[] nlocs_scat;
    //delete[] offset_scat;
    //delete[] recv_loc;
    //delete[] recv_offset;
    //delete[] TotRecvElement_IDs;
    //delete[] TotRecvElement_IDs_v;
    
    //Array<double>* Verts_loc;
    
    //std::map<int,std::vector<double> > recv_back_verts;
    
    
    return P;
}




void Example3DPartitioning(MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);

    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    int nel = 8;
    int * eltype = new int[nel];
    
    eltype[0] = 8;
    eltype[1] = 8;
    eltype[2] = 8;
    eltype[3] = 8;
    eltype[4] = 8;
    eltype[5] = 8;
    eltype[6] = 8;
    eltype[7] = 8;
    
    int npo   = 0;
    
    for(int i = 0;i < nel; i++)
    {
        npo += eltype[i];
    }
    
    int * test = new int[npo];
    
    test[0] = 0;test[1] = 1;test[2] = 6;test[3] = 5;        test[4]  = 14+0;test[5] = 14+1;test[6] = 14+6;test[7] = 14+5;
    test[8] = 1;test[9] = 2;test[10] = 7;test[11] = 6;      test[12] = 14+1;test[13] = 14+2;test[14] = 14+7;test[15] = 14+6;
    test[16] = 2;test[17] = 3;test[18] = 8;test[19] = 7;    test[20] = 14+2;test[21] = 14+3;test[22] = 14+8;test[23] = 14+7;
    test[24] = 3;test[25] = 4;test[26] = 9;test[27] = 8;    test[28] = 14+3;test[29] = 14+4;test[30] = 14+9;test[31] = 14+8;
    test[32] = 5;test[33] = 6;test[34] = 11;test[35] = 10;  test[36] = 14+5;test[37] = 14+6;test[38] = 14+11;test[39] = 14+10;
    test[40] = 6;test[41] = 7;test[42] = 12;test[43] = 11;  test[44] = 14+6;test[45] = 14+7;test[46] = 14+12;test[47] = 14+11;
    test[48] = 7;test[49] = 8;test[50] = 13;test[51] = 12;  test[52] = 14+7;test[53] = 14+8;test[54] = 14+13;test[55] = 14+12;
    test[56] = 8;test[57] = 9;test[58] = 14;test[59] = 13;  test[60] = 14+8;test[61] = 14+9;test[62] = 14+14;test[63] = 14+13;

    
    int nloc     = int(nel/world_size) + ( world_rank < nel%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(nel/world_size) + MIN(world_rank, nel%world_size);
    int* elmdist = new int[world_size];

    int npo_loc=0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += eltype[offset+i];
    }
    
    int* locs        = new int[world_size];
    int* npo_locs    = new int[world_size];
    int* npo_offset  = new int[world_size+1];
    npo_offset[0]=0;
    
    for(int i=0;i<world_size;i++)
    {
        if (i==world_rank)
        {
            locs[i]     = nloc;
            npo_locs[i] = npo_loc;
        }
        else
        {
            locs[i]     = 0;
            npo_locs[i] = 0;
        }
    }
    
    for(int i=0;i<world_size+1;i++)
    {
        if (i==world_rank)
        {
            elmdist[i]    = offset;
        }
        else
        {
            elmdist[i]    = 0;
        }
    }
    
    int* reduced_locs       = new int[world_size];
    int* reduced_npo_locs   = new int[world_size];
    int* reduced_elmdist    = new int[world_size+1];
    
    for(int i=0;i<world_size;i++)
    {
        reduced_locs[i]    = 0;
        reduced_elmdist[i] = 0;
    }
    
    MPI_Allreduce(locs,     reduced_locs,     world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(npo_locs, reduced_npo_locs, world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(elmdist,  reduced_elmdist,  world_size+1, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    
    for(int i=0;i<world_size;i++)
    {
        npo_offset[i+1] = npo_offset[i]+reduced_npo_locs[i];
    }
    
    reduced_elmdist[world_size] = nel;
    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    //std::cout << "nloc = " << nloc << std::endl;
    eptr[0]  = 0;
    
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+eltype[offset+i];
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = test[npo_offset[world_rank]+j];
            //std::cout << eind[j] << " ";
        }
        //std::cout << std::endl;
    }
    
    /*
    //map< pair<int, int>, HalfEdge* > HE = GetHalfEdges(test,eptr,nloc,offset);
    //map< pair<int, int>, HalfEdge* >::iterator it;
    

    if (world_rank == 0)
    {
        int cnt = 0;
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        for(it = HE.begin(); it != HE.end(); it++)
        {
            //std::cout << nloc << " " << cnt << std::endl;
            std::cout << cnt << " " << it->first.first << " " << it->first.second << " (" << it->second->oppositeHalfEdge->vertex << " " << it->second->oppositeHalfEdge->opposite_vertex <<")" << std::endl;
            cnt++;
            
        }
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    */
    /*
    if (world_rank == 2)
    {
        std::cout << "==========================================="<<std::endl;
        for(int i=0;i<nloc;i++)
        {
            for(int j=eptr[i];j<eptr[i+1];j++)
            {
                int vid = eind[j];
                std::cout << vid << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "==========================================="<<std::endl;
    }
    */
    /*
    if (world_rank == 3)
    {
        for(int i=0;i<nloc;i++)
        {
            for(int j=eptr[i];j<eptr[i+1];j++)
            {
                std::cout << eind[j] << " ";
            }
            
            std::cout << std::endl;
        }
    }
    */
    //===================================================================================================================================
    
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
    idx_t *adjwgt    = NULL;
    idx_t *vsize     = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt = NULL;
    
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }
    
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    
    idx_t part_[]    = {nloc};
    idx_t *part      = part_;
    
    
    ParMETIS_V3_PartMeshKway(reduced_elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
//    if(world_rank == 1
//       )
//    {
//        for(int i = 0; i < nloc; i++)
//        {
//            std::cout << part[i] << std::endl;
//        }
//    }
    
    
    ParMETIS_V3_Mesh2Dual(reduced_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
    ParMETIS_V3_AdaptiveRepart(reduced_elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
    int rank = 0;
    if(world_rank == 1)
    {
        std::cout << std::endl;
        
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    //===================================================================================================================================
    
    if(world_rank == rank)
    {
    
        //std::cout << "rank :: " << world_rank << std::endl;
        for(int i=0;i<nloc+1;i++)
        {
            std::cout << xadj[i] << " --> ";
        }
        
        for(int j=0;j<xadj[nloc];j++)
        {
            std::cout << adjncy[j] << " ";
        }
    }
}


ParVar* CreateParallelData(int N, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nloc             = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(N/size) + MIN(rank, N%size);
    
    
    int* proc_nlocs                 = new int[size];
    int* proc_offset                = new int[size];
    int* reduced_proc_nlocs             = new int[size];
    int* reduced_proc_offset            = new int[size];
    
    for(int i=0;i<size;i++)
    {
        reduced_proc_nlocs[i] = 0;
        reduced_proc_offset[i] = 0;
        
        if(i==rank)
        {
            proc_nlocs[i]  = nloc;
            proc_offset[i] = offset;
        }
        else
        {
            proc_nlocs[i]  = 0;
            proc_offset[i] = 0;
        }
    }
    
    MPI_Allreduce(proc_nlocs,  reduced_proc_nlocs,  size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(proc_offset, reduced_proc_offset, size, MPI_INT, MPI_SUM, comm);
    
    ParVar* pv = new ParVar;
    
    pv->size    = size;
    pv->nlocs   = reduced_proc_nlocs;
    pv->offsets = reduced_proc_offset;
    
    return pv;
}



// This function computes the eptr and eind array which are required for most ParMetis APIs.
// This is done based on the global element2node (e2n) array, the number of elements, the communicator
// and the type of elements we are currently using. For now the type of elements is fixed to an integer
// assuming that we use one type of element throughout the whole mesh. However this needs to become
// an int* in order to allow for hybrid meshes.
// e2n has the Nvert per element stored consecutively for each element. Hence this array is Nel*NvertPerElement long.

ParVar_ParMetis* CreateParallelDataParmetis(ParArray<int>* e2n, MPI_Comm comm, int type)
{
    int i,j;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int Nel = e2n->getNglob();
    //std::cout << "number of elements = " << Nel;
    int nloc             = int(Nel/size) + ( rank < Nel%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(Nel/size) + MIN(rank, Nel%size);
    
    int npo_loc = 0;
    for(i=0;i<nloc;i++)
    {
        npo_loc += type;
    }
    
    int* nlocs                 = new int[size];
    int* reduced_nlocs             = new int[size];
    int* npo_locs              = new int[size];
    int* reduced_npo_locs          = new int[size];

    for(i=0;i<size;i++)
    {
        reduced_nlocs[i]        = 0;
        reduced_npo_locs[i]     = 0;
        
        if(i==rank)
        {
            nlocs[i]        = nloc;
            npo_locs[i]     = npo_loc;
        }
        else
        {
            nlocs[i]        = 0;
            npo_locs[i]     = 0;
        }
    }
    
    int* elm_dist              = new int[size+1];
    int* npo_offset            = new int[size+1];
    int* reduced_elm_dist          = new int[size+1];
    int* reduced_npo_offset        = new int[size+1];
    
    for(i=0;i<size+1;i++)
    {
        reduced_elm_dist[i]   = 0;
        reduced_npo_offset[i] = 0;
        if(i==rank)
        {
            elm_dist[i]   = offset;
            npo_offset[i] = offset*type;
        }
        else
        {
            elm_dist[i]  = 0;
            npo_offset[i] = 0;
        }
    }
    
    
    MPI_Allreduce(nlocs,        reduced_nlocs,      size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs,     reduced_npo_locs,   size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(elm_dist,     reduced_elm_dist,   size+1,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_offset,   reduced_npo_offset, size+1,   MPI_INT, MPI_SUM, comm);

    reduced_elm_dist[size] = Nel;
    reduced_npo_offset[size] = Nel*type;

    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    eptr[0]  = 0;
    for(i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+type;
        
        for(j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = e2n->data[j];
        }
    }
    
    ParVar_ParMetis* pv_parmetis = new ParVar_ParMetis;
    
    pv_parmetis->size        =  size;
    pv_parmetis->nlocs       =  reduced_nlocs;
    pv_parmetis->elmdist     =  reduced_elm_dist;
    pv_parmetis->npo_locs    =  reduced_npo_locs;
    pv_parmetis->npo_offset  =  reduced_npo_offset;
    pv_parmetis->eptr        =  eptr;
    pv_parmetis->eind        =  eind;
    
//    delete[] elm_dist;
//    delete[] npo_offset;
//    delete[] reduced_elm_dist;
//    delete[] reduced_npo_offset;
    
    //delete[] eptr;
    //delete[] eind;
    
    return pv_parmetis;
}


void DivideElements(Array<int>* part_on_root, Array<int>* ien_on_root, Array<double>* xcn_on_root, MPI_Comm comm)
{
    int i=0;
    int j=0;
    int k=0;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::vector< std::vector<int> > divider(size);
    std::vector< std::vector<int> > divider_vert(size);

    int* divider_cnt   = new int[size];
    int* divider_cnt_v = new int[size];
    int* divider_cnt_c = new int[size];

    int* jag_offset  = new int[size];
    int* tmp_cntr    = new int[size];
    int* tmp_cntr2   = new int[size];
    int* tmp_cntr3   = new int[size];
    
    JagArray<int>* jArrElm = NULL;
    JagArray<int>* jArrVrt = NULL;
    JagArray<double>* jArrCrd = NULL;
    
    
    
    int Nel = part_on_root->getNrow();
    int* jag_array = NULL;
    int* jag_array_v = NULL;
    if( rank == 0)
    {
        for(i=0;i<size;i++)
        {
            divider_cnt[i]  =0;
            divider_cnt_v[i]=0;
            divider_cnt_c[i]=0;
            jag_offset[i]   =0;
            tmp_cntr[i]     =0;
            tmp_cntr2[i]    =0;
            tmp_cntr3[i]    =0;
        }
        for(i=0;i<Nel;i++)
        {
            divider_cnt[part_on_root->getVal(i,0)]=divider_cnt[part_on_root->getVal(i,0)]+1;
        }
        
        for(i=0;i<size;i++)
        {
            divider_cnt_v[i] = divider_cnt[i]*8;
            divider_cnt_c[i] = divider_cnt[i]*8*3;
        }
        
        for(i=0;i<size-1;i++)
        {
            jag_offset[i+1]=jag_offset[i]+divider_cnt[i];
        }
        
        jArrElm = new JagArray<int>(size,divider_cnt);
        jArrVrt = new JagArray<int>(size,divider_cnt_v);
        jArrCrd = new JagArray<double>(size,divider_cnt_c);
        
        int part_id = 0;
        int offset  = 0;
        int tel     = 0;
        int tel2    = 0;
        int tel3    = 0;
        jag_array   = new int[Nel];
        jag_array_v = new int[Nel*8];
        for(i=0;i<Nel;i++)
        {
            part_id = part_on_root->getVal(i,0);
            offset  = jag_offset[part_id];
            tel     = tmp_cntr[part_id];
            jag_array[offset+tel] = i;

            jArrElm->setVal(part_id,tel,i);
            for(j=0;j<8;j++)
            {
                tel2 = tmp_cntr2[part_id];
                jag_array_v[offset*8+tel2] = ien_on_root->getVal(i,j);
                jArrVrt->setVal(part_id,tel2,ien_on_root->getVal(i,j));
                tmp_cntr2[part_id] = tmp_cntr2[part_id] + 1;
                
                for(k=0;k<3;k++)
                {
                    tel3 = tmp_cntr3[part_id];
                    jArrCrd->setVal(part_id,tel3,xcn_on_root->getVal(ien_on_root->getVal(i,j),k));
                    tmp_cntr3[part_id] = tmp_cntr3[part_id] + 1;
                }
            }
            
            tmp_cntr[part_id]=tmp_cntr[part_id]+1;
        }
        
        for(i=1;i<size;i++)
        {
            /*
            int n_size = divider[i].size();
            MPI_Send(&n_size, 1, MPI_INT, i, 5678, comm);
            MPI_Send(&divider[i][0], n_size, MPI_INT, i, 5678+100, comm);
            
            int n_size_vert = divider_vert[i].size();
            MPI_Send(&n_size_vert, 1, MPI_INT, i, 5678+200, comm);
            MPI_Send(&divider_vert[i][0], n_size_vert, MPI_INT, i, 5678+300, comm);
            */
            int n_size2 = divider_cnt[i];
            MPI_Send(&n_size2, 1, MPI_INT, i, 5678+400, comm);
            MPI_Send(&jag_array[jag_offset[i]], n_size2, MPI_INT, i, 5678+500, comm);
            
            int n_size3 = divider_cnt[i]*8;
            MPI_Send(&n_size3, 1, MPI_INT, i, 5678+600, comm);
            MPI_Send(&jag_array_v[jag_offset[i]*8], n_size3, MPI_INT, i, 5678+700, comm);

            int n_size4 = divider_cnt_c[i];
            MPI_Send(&n_size4, 1, MPI_INT, i, 5678+800, comm);
            MPI_Send(&jArrCrd->data[jArrCrd->getRowOffset(i)], n_size4, MPI_INT, i, 5678+900, comm);

        }
    }
    else
    {
        /*
        int n_recv;
        MPI_Recv(&n_recv, 1, MPI_INT, 0, 5678, comm, MPI_STATUS_IGNORE);
        int* recv_el = new int[n_recv];
        MPI_Recv(&recv_el[0], n_recv, MPI_INT, 0, 5678+100, comm, MPI_STATUS_IGNORE);
        
        int n_recv_vert;
        MPI_Recv(&n_recv_vert, 1, MPI_INT, 0, 5678+200, comm, MPI_STATUS_IGNORE);
        int* recv_vert = new int[n_recv_vert];
        MPI_Recv(&recv_vert[0], n_recv_vert, MPI_INT, 0, 5678+300, comm, MPI_STATUS_IGNORE);
        */
        int n_recv2;
        MPI_Recv(&n_recv2, 1, MPI_INT, 0, 5678+400, comm, MPI_STATUS_IGNORE);
        int* recv_el2 = new int[n_recv2];
        MPI_Recv(&recv_el2[0], n_recv2, MPI_INT, 0, 5678+500, comm, MPI_STATUS_IGNORE);
        
        int n_recv3;
        MPI_Recv(&n_recv3, 1, MPI_INT, 0, 5678+600, comm, MPI_STATUS_IGNORE);
        int* recv_vert2 = new int[n_recv3];
        MPI_Recv(&recv_vert2[0], n_recv3, MPI_INT, 0, 5678+700, comm, MPI_STATUS_IGNORE);
        
        int n_recv4;
        MPI_Recv(&n_recv4, 1, MPI_INT, 0, 5678+800, comm, MPI_STATUS_IGNORE);
        int* recv_vert4 = new int[n_recv4];
        MPI_Recv(&recv_vert4[0], n_recv4, MPI_INT, 0, 5678+900, comm, MPI_STATUS_IGNORE);
        
       
        
    }
    
    delete[] jag_array;
    delete[] jag_array_v;
    
    delete jArrVrt;
    delete jArrCrd;
    delete jArrElm;
    
    delete[] divider_cnt;
    delete[] divider_cnt_c;
    delete[] divider_cnt_v;
    
    delete[] jag_offset;
    delete[] tmp_cntr;
    delete[] tmp_cntr2;
    delete[] tmp_cntr3;

    //delete divider;
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout << "rank = " << world_rank << " dividing = " << duration << std::endl;
}
