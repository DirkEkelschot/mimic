#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#include <iostream>
#include "adapt_partition.h"

using namespace std;




int FindRank(int* arr, int size, int val)
{
    int start = 0;
    int last  = size-1;
    
    int mid   = (start+last)/2;
    
    while (start<=last)
    {
        if (arr[mid]<val)
        {
            start = mid + 1;
        }
        else
        {
            last  = mid - 1;
        }
        mid = (start+last)/2;
    }
        
    return mid;
}



GathervObject* GetGathervObject(int nloc, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int* locs     = new int[size];
    int* red_locs = new int[size];

    for(int i=0;i<size;i++)
    {
        red_locs[i]  = 0;
        
        if(i==rank)
        {
            locs[i]  = nloc;
        }
        else
        {
            locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(locs, red_locs, size, MPI_INT, MPI_SUM, comm);
    
    int* red_offsets = new int[size];
    red_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        red_offsets[i+1]=red_offsets[i]+red_locs[i];
    }
    
    int length = red_offsets[size-1]+red_locs[size-1];
    
    GathervObject* gObj = new GathervObject;
    gObj->data = new int[length];
    gObj->nlocs = red_locs;
    gObj->offsets = red_offsets;
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
    int* red_exter_nlocs             = new int[size];
    int* red_exter_offset            = new int[size];
    

    
    for(int i=0;i<size;i++)
    {
        red_exter_nlocs[i]  = 0;
        red_exter_offset[i] = 0;
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
                  red_exter_nlocs,
                  size,
                  MPI_INT,
                  MPI_SUM,
                  comm);
    
    
    red_exter_offset[0]=0;
        
    
    //
    for(int i=0;i<size-1;i++)
    {
        red_exter_offset[i+1]=red_exter_offset[i]+red_exter_nlocs[i];
    }
    //
    
    int nexter_tot = red_exter_offset[size-1]+red_exter_nlocs[size-1];
    
    /*
    std::vector<int> recv(nexter_tot);
     
    MPI_Allgatherv(&u_faces[0],
                   exter_size,
                   MPI_INT,
                   &recv[0],
                   red_exter_nlocs,
                   red_exter_offset,
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
    delete[] red_exter_nlocs;
    delete[] red_exter_offset;
    
    
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
    int* red_locs = new int[size];

    for(int i=0;i<size;i++)
    {
        red_locs[i]  = 0;
        
        if(i==rank)
        {
            locs[i]  = locvec.size();
        }
        else
        {
            locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(locs, red_locs, size, MPI_INT, MPI_SUM, comm);
    
    int* red_offsets = new int[size];
    red_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        red_offsets[i+1]=red_offsets[i]+red_locs[i];
    }
    
    int tot = red_offsets[size-1]+red_locs[size-1];
    
    ParArrayOnRoot* parr_root = new ParArrayOnRoot;
    parr_root->data           = new int[tot];
    parr_root->nlocs          = red_locs;
    parr_root->offsets        = red_offsets;
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

Partition* CollectVerticesPerRank(ParArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm)
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
    ParArrayOnRoot* gathered_on_root = GatherVecToRoot(loc_elems, comm);
    
    
    int* nlocs = new int[size];
    int* offset = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs[i] = gathered_on_root->nlocs[i]*3;
        offset[i] = gathered_on_root->offsets[i]*3;
    }
    
    Partition* parti = new Partition;
    
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
        verts = new double[gathered_on_root->length*3];

        for(int i=0;i<gathered_on_root->length;i++)
        {
            verts[i*3+0] = xcn_r->getVal(gathered_on_root->data[i],0);
            verts[i*3+1] = xcn_r->getVal(gathered_on_root->data[i],1);
            verts[i*3+2] = xcn_r->getVal(gathered_on_root->data[i],2);
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


ParArray<int>* DeterminePartitionLayout(ParArray<int>* ien, MPI_Comm comm)
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
    
    ParallelState_Parmetis* pstate_parmetis = new ParallelState_Parmetis(ien,comm,8);
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
        
    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj,&adjncy,&comm);
    
    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj,
                         adjncy,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part, &comm);
     
    ParArray<int>*  part_arr = new ParArray<int>(ien->getNglob(),1,comm);
    part_arr->data = part;
    
    delete pstate_parmetis;
    
    return part_arr;
    
}


// This function determines a map that gets the unique list of elements for that need to be requested from a given rank other than current rank.
Array<double>* DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* part, Array<double>* xcn_on_root, ParArray<double>* xcn, MPI_Comm comm)
{
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
        
    int el_id;
    int p_id;
    int v_id; 
    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > verts_to_send_to_ranks;   

    set<int> u_verts_set;
    std::vector<int> u_verts_vec;

    std::map<int,set<int> > u_verts_other_ranks_set;
    std::map<int,std::vector<int> > u_verts_other_ranks_vec;
    std::map<int,int> cnt_other_ranks;
    
    std::vector<int> part_v;
    int r = 0;
    int l_id=0;
    for(i=0;i<part->getNrow();i++)
    {
        p_id  = part->getVal(i,0);
        el_id = part->getParallelState()->getOffset(rank)+i;

        if(p_id!=rank)
        {
            elms_to_send_to_ranks[p_id].push_back(el_id);
            
            for(int k=0;k<8;k++)
            {
                v_id = ien->getVal(i,k);
                verts_to_send_to_ranks[p_id].push_back(v_id);
            }
        }
        else
        {
            for(int k=0;k<8;k++)
            {
                v_id = ien->getVal(i,k);
                
                if(u_verts_set.find( v_id ) == u_verts_set.end())
                {
                    u_verts_set.insert(v_id);
                    u_verts_vec.push_back(v_id);
                    
                    r = FindRank(xcn->getParallelState()->getOffsets(),size+1,v_id);
                    part_v.push_back(r);
                    
                    l_id++;
                }
            }
        }
    }
    
    int to_send_size = elms_to_send_to_ranks.size();
    
    int* red_to_send_size = new int[size];
    int* arr_to_send_size = new int[size];
    
    for(i=0;i<size;i++)
    {
        red_to_send_size[i] = 0;
        
        if(i==rank)
        {
            arr_to_send_size[i] = to_send_size+1;
        }
        else
        {
            arr_to_send_size[i] = 0;
        }
    }

    MPI_Allreduce(arr_to_send_size, red_to_send_size, size, MPI_INT, MPI_SUM, comm);
    
    int* red_to_send_size_offset = new int[size];
    int offset = 0;
    for(i=0;i<size;i++)
    {
        red_to_send_size_offset[i] = offset;
        offset = offset+red_to_send_size[i];
    }
    int send_map_size = 0;
    int loc_size = to_send_size+1; // This size is added by one since we add the rank number to the array.
    
    MPI_Allreduce(&loc_size, &send_map_size, 1, MPI_INT, MPI_SUM, comm);
    
    int* send_map_part_id            = new int[send_map_size];
    int* send_map_Nel_per_part_id    = new int[send_map_size];
    for(i=0;i<send_map_size;i++)
    {
        send_map_part_id[i]         = 0;
        send_map_Nel_per_part_id[i] = 0;
    }
    int* from_to_rank       = new int[to_send_size+1];
    int* num_elms_to_rank   = new int[to_send_size+1];
    from_to_rank[0]         = rank;
    num_elms_to_rank[0]     = -1;
    int t = 1;
    
    std::map<int,std::vector<int> >::iterator it;
    for(it=elms_to_send_to_ranks.begin();it!=elms_to_send_to_ranks.end();it++)
    {
        from_to_rank[t]     = it->first;
        num_elms_to_rank[t] = it->second.size();
        t++;
    }

    MPI_Allgatherv(&from_to_rank[0], loc_size, MPI_INT,
                   &send_map_part_id[0],
                   red_to_send_size, red_to_send_size_offset, MPI_INT,comm);
    
    MPI_Allgatherv(&num_elms_to_rank[0], loc_size, MPI_INT,
                   &send_map_Nel_per_part_id[0],
                   red_to_send_size, red_to_send_size_offset, MPI_INT,comm);
    
    std::map<int, set<int> > s_iset_map;
    std::map<int, std::vector<int> > r_ivec_map;
    std::map<int, std::vector<int> > r_ivec_size_map;
    
    //=========================================
    for(i=0;i<size;i++)
    {
        int of = red_to_send_size_offset[i];
        int nl = red_to_send_size[i];
        for(int j=of+1;j<of+nl;j++)
        {
            s_iset_map[send_map_part_id[of]].insert(send_map_part_id[j]);
            r_ivec_map[send_map_part_id[j]].push_back(send_map_part_id[of]);
            r_ivec_size_map[send_map_part_id[j]].push_back(send_map_Nel_per_part_id[j]);
        }
    }

    std::map<int,std::map<int,int> > loc_alloc;
    
    for(it=r_ivec_map.begin();it!=r_ivec_map.end();it++)
    {
        for(int k=0;k<it->second.size();k++)
        {
            loc_alloc[it->first][r_ivec_map[it->first][k]] = r_ivec_size_map[it->first][k];
        }
    }
    
    //============================================================================
    //================ Print the receiving schedule for now;======================
    //============================================================================
//    if(rank == 0)
//    {
//        for(it=r_ivec_map.begin();it!=r_ivec_map.end();it++)
//        {
//            std::cout << "rank " << it->first << " receives an array of size ";
//            for(int k=0;k<it->second.size();k++)
//            {
//                std::cout << r_ivec_size_map[it->first][k] << " from " << r_ivec_map[it->first][k] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    //============================================================================
    //================ Print the receiving schedule for now;======================
    //============================================================================
    
    
    std::map<int,int> alloc = loc_alloc[rank];
    
    int* recv_offset = new int[loc_alloc.size()+1];
    recv_offset[0]   = 0;
    int* recv_loc    = new int[loc_alloc.size()];;
    int recv_size    = 0;
    
    std::map< int, int> offset_map;
    std::map< int, int> loc_map;
    std::map< int, int>::iterator it_loc;
    i = 0;
    for(it_loc=alloc.begin();it_loc!=alloc.end();it_loc++)
    {
        recv_loc[i]      = it_loc->second;
        recv_offset[i+1] = recv_offset[i]+recv_loc[i];
        recv_size = recv_size+it_loc->second;

        loc_map[it_loc->first]=recv_loc[i];
        offset_map[it_loc->first]=recv_offset[i];

        i++;
    }
    int* recv_collector   = new int[recv_size];
    int* recv_collector_v = new int[recv_size*8];
    
    
    for(int i =0;i<recv_size;i++)
    {
        recv_collector[i] = 10;
    }
    
    std::map< int, std::map< int, int> > s_recv_alloc;
    
    int n_req_recv;
    
    int n_req_recv_v;
    for(int q=0;q<size;q++)
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
                
                MPI_Send(&n_req_v, 1, MPI_INT, dest, 9000+dest, comm);
                MPI_Send(&verts_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                
                i++;
            }
        }
        else if (s_iset_map[q].find( rank ) != s_iset_map[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_collector[offset_map[q]], n_req_recv, MPI_INT, q, 100+rank*2, comm, MPI_STATUS_IGNORE);
            
            
            
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, 9000+rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_collector_v[offset_map[q]*8], n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
        }
    }
    
    
    
    
    
    
    std::vector<int> u_new_verts_vec;
    std::map<int,std::vector<int> > crds_to_send_to_ranks;

    int v = 0;
    for(int i=0;i<recv_size*8;i++)
    {
        int v_id_n = recv_collector_v[i];
        if(u_verts_set.find( v_id_n ) == u_verts_set.end())
        {
            u_verts_set.insert(v_id_n);
            u_verts_vec.push_back(v_id_n);
            
            r = FindRank(xcn->getParallelState()->getOffsets(),size,v_id_n);
            part_v.push_back(r);
            
            if (r!=rank)
            {
                crds_to_send_to_ranks[r].push_back(v_id_n);
                //std::cout << v_id_n << " " << r << std::endl;
                
            }
            v++;
        }
    }
    
    std::map<int,std::vector<int> >::iterator it4;
    for(it4=crds_to_send_to_ranks.begin();it4!=crds_to_send_to_ranks.end();it4++)
    {
        std::cout << rank<< " send " << it4->second.size() << " array with values: ";
        std::vector<int>::iterator it5;
        for(it5=it4->second.begin();it5!=it4->second.end();it5++)
        {
            std::cout << *it5 << " ";
        }
        
        std::cout << "to rank " << it4->first << std::endl;
    }
    
    int crds_to_send_size = crds_to_send_to_ranks.size();
    
    int* red_crds_to_send_size = new int[size];
    int* arr_crds_to_send_size = new int[size];
    
    for(i=0;i<size;i++)
    {
        red_crds_to_send_size[i] = 0;
        
        if(i==rank)
        {
            arr_crds_to_send_size[i] = crds_to_send_size+1;
        }
        else
        {
            arr_crds_to_send_size[i] = 0;
        }
    }

    MPI_Allreduce(arr_crds_to_send_size, red_crds_to_send_size, size, MPI_INT, MPI_SUM, comm);
    
    int* red_crds_to_send_size_offset = new int[size];
    offset = 0;
    for(i=0;i<size;i++)
    {
        red_crds_to_send_size_offset[i] = offset;
        offset = offset+red_crds_to_send_size[i];
    }
    int crds_send_map_size = 0;
    int crds_loc_size = crds_to_send_size+1; // This size is added by one since we add the rank number to the array.
    
    MPI_Allreduce(&crds_loc_size, &crds_send_map_size, 1, MPI_INT, MPI_SUM, comm);
    //std::cout << "crds_send_map_size " << crds_send_map_size << std::endl;
    int* crds_send_map_part_id            = new int[crds_send_map_size];
    int* send_map_Nverts_per_part_id      = new int[crds_send_map_size];
    for(i=0;i<crds_send_map_size;i++)
    {
        crds_send_map_part_id[i]         = 0;
        send_map_Nverts_per_part_id[i]   = 0;
    }
    int* crds_from_to_rank       = new int[crds_to_send_size+1];
    int* num_crds_to_rank        = new int[crds_to_send_size+1];
    crds_from_to_rank[0]         = rank;
    num_crds_to_rank[0]          = -1;
    t = 1;
    
    std::map<int,std::vector<int> >::iterator it2;
    for(it2=crds_to_send_to_ranks.begin();it2!=crds_to_send_to_ranks.end();it2++)
    {
        crds_from_to_rank[t]     = it2->first;
        num_crds_to_rank[t]      = it2->second.size();
        t++;
    }

    MPI_Allgatherv(&crds_from_to_rank[0], crds_loc_size, MPI_INT,
                   &crds_send_map_part_id[0],
                   red_crds_to_send_size, red_crds_to_send_size_offset, MPI_INT,comm);
    
    MPI_Allgatherv(&num_crds_to_rank[0], crds_loc_size, MPI_INT,
                   &send_map_Nverts_per_part_id[0],
                   red_crds_to_send_size, red_crds_to_send_size_offset, MPI_INT,comm);
    
//    if(rank == 1)
//    {
//        for(i=0;i<crds_send_map_size;i++)
//        {
//            std::cout << i << " " << crds_send_map_part_id[i] << " " << send_map_Nverts_per_part_id[i] << std::endl;
//        }
//        //std::cout << std::endl;
//    }
    
//    std::map<int, set<int> > s_iset_map;
//    std::map<int, std::vector<int> > r_ivec_map;
//    std::map<int, std::vector<int> > r_ivec_size_map;
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    // ---> First attempt <--- //
    
    ParArrayOnRoot* gathered_on_root = GatherVecToRoot(u_verts_vec, comm);
    
    int* nlocs_scat = new int[size];
    int* offset_scat = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs_scat[i]  = gathered_on_root->nlocs[i]*3;
        offset_scat[i] = gathered_on_root->offsets[i]*3;
    }
       
    double* verts = NULL;
    Array<double>* Verts_loc = new Array<double>(u_verts_vec.size(),3);
    if(rank == 0)
    {
        
        verts = new double[gathered_on_root->length*3];
        for(int i=0;i<gathered_on_root->length;i++)
        {
            verts[i*3+0] = xcn_on_root->getVal(gathered_on_root->data[i],0);
            verts[i*3+1] = xcn_on_root->getVal(gathered_on_root->data[i],1);
            verts[i*3+2] = xcn_on_root->getVal(gathered_on_root->data[i],2);
        }
        
    }
    
    MPI_Scatterv(&verts[0], nlocs_scat, offset_scat, MPI_DOUBLE, &Verts_loc->data[0], u_verts_vec.size()*3, MPI_DOUBLE, 0, comm);
    
    
    delete[] verts;
    delete[] nlocs_scat;
    delete[] offset_scat;
    delete[] recv_loc;
    delete[] recv_offset;
    delete[] recv_collector;
    delete[] recv_collector_v;
    
    
    
    
    return Verts_loc;
}


Partition* CollectElementsPerRank(ParArray<int>* ien, Array<int>* ien_root, MPI_Comm comm)
{
    int i;
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
    
    ParallelState_Parmetis* pstate_parmetis = new ParallelState_Parmetis(ien,comm,8);
//
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
    
    //idx_t *elmwgt;

//    int np           = size;
//    idx_t ncon_[]    = {1};
//    idx_t *ncon      = ncon_;
//    real_t *tpwgts   = new real_t[np*ncon[0]];
//
//    for(i=0; i<np*ncon[0]; i++)
//    {
//        tpwgts[i] = 1.0/np;
//    }
    int* part = new int[nloc];
    
    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj,&adjncy,&comm);
    /*
    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj,
                         adjncy,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                        &edgecut, part, &comm);
    
    */
    
    
    /* 
    for(int i=0;i<nloc;i++)
    {
	std::cout << rank << " :: " << i+pstate_parmetis->getElmdistAtRank(rank) << " ---> "; 
	for(int j=xadj[i];j<xadj[i+1];j++)
	{
		std::cout << adjncy[j] << " ";
	}
	std::cout << std::endl;
    }
    
    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(), 
			      xadj, adjncy, 
			      elmwgt, adjwgt, vsize, 
			      wgtflag, numflag, ncon, 
			      nparts, tpwgts, ubvec, 
			      itr, options, &edgecut, part, &comm);
    */
    /*
    ParMETIS_V3_PartMeshKway(pstate_parmetis->getElmdist(),
                             pstate_parmetis->getEptr(),
                             pstate_parmetis->getEind(),
                             elmwgt, wgtflag, numflag,
                             ncon, ncommonnodes, nparts,
                             tpwgts, ubvec, options,
                             &edgecut, part, &comm);
    */
    /*
    
    if(rank == 0 || rank == 1)
    {
	for(int i=0;i<ien->getNrow();i++)
	{       std::cout <<rank <<  " Element # = " << i+pstate_parmetis->getElmdistAtRank(rank) << " :: ";
		for(int j=0;j<ien->getNcol();j++)
		{
			std::cout << ien->getVal(i,j) << " ";
		}
		std::cout << std::endl;
	}
    }
    */
    //std::cout << rank << " nloc " << nloc << std::endl; 
    //if(rank == 1)
    //{
    int cnt = 0;
	for(i=0;i<nloc;i++)
	{
		if(part[i]!=rank)
	        {
		    //std::cout << "rank  " << rank << " should send " << i+pstate_parmetis->getElmdistAtRank(rank) << " to rank " << part[i]  << std::endl;
			cnt++;
		}
	}
    //}
    std::cout << rank << " needs to send " << cnt << " elements " << std::endl;
    ParArray<int>*  part_arr = new ParArray<int>(ien->getNglob(),1,comm);
    part_arr->data = part;
    int tot = ien->getNglob();
    
    Array<int>* output = NULL;
    
    if (rank == 0)
    {
        output = new Array<int>(tot,1);
    }

    MPI_Gatherv(&part_arr->data[0],
                   nloc,
                   MPI_INT,
                   &output->data[0],
                   part_arr->getParallelState()->getNlocs(),
                   part_arr->getParallelState()->getOffsets(),
                   MPI_INT, 0, comm);
    /*
    MPI_Allgatherv(&part_arr->data[0],
                     nloc,
                     MPI_INT,
                     &output->data[0],
                     part_arr->getParallelState()->getNlocs(),
                     part_arr->getParallelState()->getOffsets(),
                     MPI_INT, comm);
    */
    //ParArrayOnRoot* gRoot = GatherArrToRoot(part, nloc, comm);
    //=================================================================
    //=================================================================
    //=================================================================
    
    std::map<int,int> loc2glob;
    std::map<int,int> glob2loc;
    
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
    //Partition* parti = new Partition;
    
    ParArrayOnRoot* gathered_on_root = GatherVecToRoot(loc_elems, comm);
    
    int* nlocs  = new int[size];
    int* offset = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs[i]  = gathered_on_root->nlocs[i]*8;
        offset[i] = gathered_on_root->offsets[i]*8;
    }

    Partition* parti = new Partition;

    parti->Verts = new Array<double>(loc_elems.size(),8);

    parti->loc2glob_Vmap = loc2glob;
    parti->glob2loc_Vmap = glob2loc;

    //parti->loc2glob_Varr = loc2glob_vert;
    //parti->glob2loc_Varr = glob2loc_vert;

    parti->xadj          = xadj;
    parti->adjncy        = adjncy;
    parti->ien           = ien;
    parti->ndim          = 3;
    
    double * verts = NULL;
    
    if(rank == 0)
    {
        verts = new double[gathered_on_root->length*8];

        for(int i=0;i<gathered_on_root->length;i++)
        {
            verts[i*8+0] = ien_root->getVal(gathered_on_root->data[i],0);
            verts[i*8+1] = ien_root->getVal(gathered_on_root->data[i],1);
            verts[i*8+2] = ien_root->getVal(gathered_on_root->data[i],2);
            verts[i*8+3] = ien_root->getVal(gathered_on_root->data[i],3);
            verts[i*8+4] = ien_root->getVal(gathered_on_root->data[i],4);
            verts[i*8+5] = ien_root->getVal(gathered_on_root->data[i],5);
            verts[i*8+6] = ien_root->getVal(gathered_on_root->data[i],6);
            verts[i*8+7] = ien_root->getVal(gathered_on_root->data[i],7);
        }
    }
   
    MPI_Scatterv(&verts[0], nlocs, offset, MPI_DOUBLE, &parti->Verts->data[0], loc_elems.size()*8, MPI_DOUBLE, 0, comm);
    
    delete part_arr;
    //delete[] elmwgt;
    delete[] xadj;
    delete[] adjncy;
    delete[] verts;
    delete pstate_parmetis;
    delete gathered_on_root;
    
    return parti;
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
    
    int* red_locs       = new int[world_size];
    int* red_npo_locs   = new int[world_size];
    int* red_elmdist    = new int[world_size+1];
    
    for(int i=0;i<world_size;i++)
    {
        red_locs[i]    = 0;
        red_elmdist[i] = 0;
    }
    
    MPI_Allreduce(locs,     red_locs,     world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(npo_locs, red_npo_locs, world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(elmdist,  red_elmdist,  world_size+1, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    
    for(int i=0;i<world_size;i++)
    {
        npo_offset[i+1] = npo_offset[i]+red_npo_locs[i];
    }
    
    red_elmdist[world_size] = nel;
    
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
    
    
    ParMETIS_V3_PartMeshKway(red_elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
//    if(world_rank == 1
//       )
//    {
//        for(int i = 0; i < nloc; i++)
//        {
//            std::cout << part[i] << std::endl;
//        }
//    }
    
    
    ParMETIS_V3_Mesh2Dual(red_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
    ParMETIS_V3_AdaptiveRepart(red_elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
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
    int* red_proc_nlocs             = new int[size];
    int* red_proc_offset            = new int[size];
    
    for(int i=0;i<size;i++)
    {
        red_proc_nlocs[i] = 0;
        red_proc_offset[i] = 0;
        
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
    
    MPI_Allreduce(proc_nlocs,  red_proc_nlocs,  size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(proc_offset, red_proc_offset, size, MPI_INT, MPI_SUM, comm);
    
    ParVar* pv = new ParVar;
    
    pv->size    = size;
    pv->nlocs   = red_proc_nlocs;
    pv->offsets = red_proc_offset;
    
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
    int* red_nlocs             = new int[size];
    int* npo_locs              = new int[size];
    int* red_npo_locs          = new int[size];

    for(i=0;i<size;i++)
    {
        red_nlocs[i]        = 0;
        red_npo_locs[i]     = 0;
        
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
    int* red_elm_dist          = new int[size+1];
    int* red_npo_offset        = new int[size+1];
    
    for(i=0;i<size+1;i++)
    {
        red_elm_dist[i]   = 0;
        red_npo_offset[i] = 0;
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
    
    
    MPI_Allreduce(nlocs,        red_nlocs,      size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs,     red_npo_locs,   size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(elm_dist,     red_elm_dist,   size+1,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_offset,   red_npo_offset, size+1,   MPI_INT, MPI_SUM, comm);

    red_elm_dist[size] = Nel;
    red_npo_offset[size] = Nel*type;

    
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
    pv_parmetis->nlocs       =  red_nlocs;
    pv_parmetis->elmdist     =  red_elm_dist;
    pv_parmetis->npo_locs    =  red_npo_locs;
    pv_parmetis->npo_offset  =  red_npo_offset;
    pv_parmetis->eptr        =  eptr;
    pv_parmetis->eind        =  eind;
    
//    delete[] elm_dist;
//    delete[] npo_offset;
//    delete[] red_elm_dist;
//    delete[] red_npo_offset;
    
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
