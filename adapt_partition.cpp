#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#include <iostream>
#include "adapt_partition.h"

using namespace std;



std::vector<int> GetAdjacencyForUS3D_V4(ParallelArray<int>* ief, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    

//    Array<int>*type;
    std::vector<int> ief_copy(ief->nloc*(ief->ncol-1));
    set<int> unique_faces;
    set<int> double_faces;
    std::vector<int> d_faces;
    std::vector<int> u_faces;
    
    
    int fid;
    
    for(int i=0;i<ief->nloc;i++)
    {
        for(int j=0;j<ief->ncol-1;j++)
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
    
    //std::cout << "sizes = = " << unique_faces.size() << " " << double_faces.size() << std::endl;
    
    //std::vector<int> double_faces2 = FindDuplicates(d_faces);
    
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
    
    std::clock_t start;
    double duration;
    //start = std::clock();
    
    // copy the vector<int> into an int* for the FindDuplicatesInParallel routine.
    int* uf_arr = new int[u_faces.size()];
    for(int i=0;i<u_faces.size();i++)
    {
        uf_arr[i] = u_faces[i];
    }
    
    std::vector<int> recv2 = FindDuplicatesInParallel(uf_arr, u_faces.size(), nexter_tot, comm);
    
    /*
    if (rank == 0)
    {
        for(int i=0;i<recv2.size();i++)
        {
            std::cout << rank << " :: " << recv2[i] << std::endl;
        }
    }
    */
    /*
    int inter_size = d_faces.size();

    int* inter_nlocs                 = new int[size];
    int* inter_offset                = new int[size];
    int* red_inter_nlocs             = new int[size];
    int* red_inter_offset            = new int[size];

    for(int i=0;i<size;i++)
    {
        red_inter_nlocs[i]  = 0;
        red_inter_offset[i] = 0;
        if(i==rank)
        {
            inter_nlocs[i] = inter_size;
        }
        else
        {
            inter_nlocs[i]  = 0;
            inter_offset[i] = 0;
        }
    }
    
    MPI_Allreduce(inter_nlocs,
                  red_inter_nlocs,
                  size,
                  MPI_INT,
                  MPI_SUM,
                  comm);
    
    
    red_inter_offset[0]=0;
        
    //
    for(int i=0;i<size-1;i++)
    {
        red_inter_offset[i+1]=red_inter_offset[i]+red_inter_nlocs[i];
    }
    //
    
    int ninter_tot      = red_inter_offset[size-1]+red_inter_nlocs[size-1];
    
    int* internal_faces = new int[ninter_tot];
    
    MPI_Allgatherv(&d_faces[0],
                   d_faces.size(),
                   MPI_INT,
                   &internal_faces[0],
                   red_inter_nlocs,
                   red_inter_offset,
                   MPI_INT, comm);
    
    std::cout << "internal and total faces " << ninter_tot+recv2.size() << " " << ief->nglob*6 << " " << ief->nglob*6-ninter_tot+recv2.size() << std::endl;
    
    */
    
    return recv2;
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout << rank << " Find array = " << duration << std::endl;
    
//
//    start = std::clock();
//
//    std::vector<int> recv23 = FindDuplicatesInParallel_Vec(u_faces, u_faces.size(), nexter_tot, comm);
//
//    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//    std::cout << rank << " Find vector = " << duration << std::endl;
//
    
    //if(rank == 0)
    //{
    //    std::cout << rank << " recv2.size() V4 = " << recv2.size() << " " << std::endl;
    //    for(int i=0;i<recv2.size();i++)
    //    {
    //        std::cout << i << "  " << recv2[i] << " " << recv2.size() <<  std::endl;
    //    }
    //}
    
    
    
//std::cout << "ada " << recv2.size() << std::endl;
        
//    for(int i=0;i<recv2.size();i++)
//    {
//        std::cout << "there we go! " << rank << " " << recv2[i] << std::endl;
//    }
//
//    std::cout << recv
//
//    std::cout << "sizew = " << recv2.size() << " " << recv.size() << " " << d_faces.size() << std::endl;
//
    /*
        
        
    std::cout << "hello = " << red_exter_offset[rank] << " " << red_exter_nlocs[rank] << std::endl;
    for(int i=red_exter_offset[rank];i<red_exter_offset[rank]+red_exter_nlocs[rank];i++)
    {
        //std::cout << rank << " " << i << " " << recv[i] << std::endl;
    }
        

     
     */
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
    parr_root->data = new int[tot];
    parr_root->nlocs = red_locs;
    parr_root->offsets = red_offsets;
    parr_root->size    = size;
    parr_root->length  = tot;
    
    MPI_Gatherv(&locvec[0],
                locvec.size(),
                MPI_INT,
                &parr_root->data[0],
                parr_root->nlocs,
                parr_root->offsets,
                MPI_INT,0, comm);
    
    return parr_root;
}



Partition* CollectVerticesPerRank(ParallelArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,int> loc2glob;
    std::map<int,int> glob2loc;
    Array<int>* loc2glob_vert = new Array<int>(ien->nloc,ien->ncol-1);
    Array<int>* glob2loc_vert = new Array<int>(ien->nloc,ien->ncol-1);
    std::vector<int> loc_verts;
    set<int> loc_verts_set;
    int gid = 0;
    int nverts_req;
    int lid = 0;
    for(int i=0;i<ien->nloc;i++)
    {
        for(int j=0;j<ien->ncol-1;j++)
        {
            gid = ien->getVal(i,j+1)-1;
            loc2glob_vert->setVal(i,j,lid);
            glob2loc_vert->setVal(i,j,gid);
            if ( loc_verts_set.find( gid ) == loc_verts_set.end() )
            {
                loc_verts.push_back(gid);
                loc_verts_set.insert(gid);
                loc2glob[lid] = gid;
                glob2loc[gid] = lid;
                lid++;
            }
        }
    }
    
    // This vector is empty on all other procs except root;
    ParArrayOnRoot* gathered_on_root = GatherVecToRoot(loc_verts, comm);
    double* loc_verts2 = new double[loc_verts.size()*3];
    double * verts;
    
    int* nlocs = new int[size];
    int* offset = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs[i] = gathered_on_root->nlocs[i]*3;
        offset[i] = gathered_on_root->offsets[i]*3;
        
    }
    
    Partition* part = new Partition;
    
    part->Verts = new Array<double>(loc_verts.size(),3);
    part->loc2glob_Vmap = loc2glob;
    part->glob2loc_Vmap = glob2loc;
    
    part->loc2glob_Varr = loc2glob_vert;
    part->glob2loc_Varr = glob2loc_vert;
    
    part->ien = ien;
    part->ndim = 3;
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

    MPI_Scatterv(&verts[0], nlocs, offset, MPI_DOUBLE, &part->Verts->data[0], loc_verts.size()*3, MPI_DOUBLE, 0, comm);
    
    return part;
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
    std::cout << "nloc = " << nloc << std::endl;
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

    idx_t *elmwgt;
    
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
    
    if(world_rank == 1
       )
    {
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    
    
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

ParVar_ParMetis* CreateParallelDataParmetis(Array<int>* e2n, MPI_Comm comm, int type)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int Nel = e2n->nloc;
    //std::cout << "number of elements = " << Nel;
    int nloc             = int(Nel/size) + ( rank < Nel%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(Nel/size) + MIN(rank, Nel%size);
    
    int npo_loc = 0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += type;
    }
    
    int* nlocs                 = new int[size];
    int* red_nlocs             = new int[size];
    int* npo_locs              = new int[size];
    int* red_npo_locs          = new int[size];

    for(int i=0;i<size;i++)
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
    
    for(int i=0;i<size+1;i++)
    {
        red_elm_dist[i]   = 0;
        red_npo_offset[i] = 0;
        if(i==rank)
        {
            elm_dist[i]   = offset;
            npo_offset[i] = offset*8;
        }
        else
        {
            elm_dist[i]  = 0;
            npo_offset[i] = 0;
        }
    }
    
    
    MPI_Allreduce(nlocs,        red_nlocs,      size,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs,     red_npo_locs,   size,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(elm_dist,     red_elm_dist,   size+1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_offset,   red_npo_offset, size+1,   MPI_INT, MPI_SUM, comm);

    red_elm_dist[size] = Nel;
    red_npo_offset[size] = Nel*8;

    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    eptr[0]  = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+type;
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = e2n->data[npo_offset[rank]+j];
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
    
    return pv_parmetis;
}
