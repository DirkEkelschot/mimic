#include "adapt_partition.h"

Partition::Partition(ParArray<int>* ien, ParArray<int>* iee, ParArray<int>* ief, ParArray<int>* ie_Nv, ParArray<int>* ie_Nf, ParArray<int>* ifn, ParArray<int>* ife, ParArray<int>* if_ref, ParArray<int>* if_Nv,  ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, ParallelState* ife_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    ien_pstate = ien_parstate;
    xcn_pstate = xcn_parstate;
    ife_pstate = ife_parstate;
    // This function computes the xadj and adjcny array and the part array which determines which element at current rank should be sent to other ranks.
    NelGlob = ien->getNglob();
    double t0 = MPI_Wtime();
    // This routine essentially determines based on the current element layout what the ideal layout should be.
    DeterminePartitionLayout(ien, pstate_parmetis, comm);

    double t1 = MPI_Wtime();
    double time_layout = t1-t0;
    double max_time_layout = 0.0;
    MPI_Allreduce(&time_layout, &max_time_layout, 1, MPI_DOUBLE, MPI_MAX, comm);

    eloc = 0;
    vloc = 0;
    floc = 0;

    // This function takes care of the send and receive operations in order to send the appropriate elements and corresponding vertices to the appropriate rank.
    // These operations are based on the fact that part holds the desired partitioning of the elements. the spread of the vertices is based on the fact that all the vertices stored in xcn are distributed "uniformly";

    //Sends each element and its vertices and faces to the correct proc.
    DetermineElement2ProcMap(ien, ief, ie_Nv, ie_Nf, xcn, U, comm);
//
    iee_part_map = getElement2EntityPerPartition(iee,  Loc_Elem_Nf,   comm);
    ief_part_map = getElement2EntityPerPartition(ief,  Loc_Elem_Nf,   comm);
    ien_part_map = getElement2EntityPerPartition(ien,  Loc_Elem_Nv,   comm);
//
    if_Nv_part_map      = getFace2EntityPerPartition(if_Nv ,   comm);
    ifn_part_map        = getFace2NodePerPartition(ifn     ,   comm);

    //ifn_part_map        = getFace2EntityPerPartition(ifn   ,   comm);
    ife_part_map        = getFace2EntityPerPartition(ife   ,   comm);
    if_ref_part_map     = getFace2EntityPerPartition(if_ref,   comm);
//
    DetermineAdjacentElement2ProcMapUS3D(ien, iee_part_map->i_map, part, xcn, U, comm);
//
    CreatePartitionDomain();
//
    nLocAndAdj_Elem = LocAndAdj_Elem.size();
 
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void Partition::DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, MPI_Comm comm)
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
    idx_t ncommonnodes_[] = {pstate_parmetis->getNcommonNodes()};
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
    
    
//    for(int u=0;u<nloc;u++)
//    {
//        part_arr[u] = rank;
//    }
    

//    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(),
//                           xadj_par, adjncy_par,
//                               elmwgt, adjwgt,
//                       vsize, wgtflag,
//                   numflag, ncon, nparts,
//                   tpwgts, ubvec, itr, options,
//                   &edgecut, part_arr, &comm);
    
    
    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);
    
    
    part = new ParArray<int>(ien->getNglob(),1,comm);
    part_global = new Array<int>(ien->getNglob(),1);

    part->data = part_arr;
//    xadj = xadj_par;
//    adjcny = adjncy_par;
    
    MPI_Allgatherv(&part->data[0],
                   nloc, MPI_INT,
                   &part_global->data[0],
                   ien_pstate->getNlocs(),
                   ien_pstate->getOffsets(),
                   MPI_INT,comm);
    
    
    delete[] xadj_par;
    delete[] adjncy_par;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void Partition::DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParArray<int>* ie_Nv, ParArray<int>* ie_Nf, ParArray<double>* xcn, Array<double>* U, MPI_Comm comm)
{
    int floc_tmp=0;
    int vloc_tmp=0;
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
    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceIDs_to_send_to_ranks;
    std::map<int,std::vector<double> > varia_to_send_to_ranks;
    

    
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_face;
    std::vector<int> faceIDs_on_rank;
    

    
        
    std::vector<int> vertIDs_on_rank;
    std::vector<int> part_v;
    
    

    
    int r     = 0;
    int lv_id = 0;
    int lf_id = 0;
    int f_id  = 0;


    int xcn_o = xcn->getOffset(rank);
    
    int ien_o = part->getOffset(rank);
    double varia = 0.0;
    int not_on_rank=0;
    int on_rank = 0;
    int* new_V_offsets = new int[size];
    int* new_F_offsets = new int[size];
    for(i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
        new_F_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    int nvPerEl;
    int nfPerEl;
    for(i=0;i<part->getNrow();i++)
    {
        p_id    = part->getVal(i,0);
        el_id   = part->getOffset(rank)+i;
        varia   = U->getVal(i,0);
        nvPerEl = ie_Nv->getVal(i,0);
        nfPerEl = ie_Nf->getVal(i,0);
        
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            nvPerElms_to_send_to_ranks[p_id].push_back(nvPerEl);
            nfPerElms_to_send_to_ranks[p_id].push_back(nfPerEl);

            varia_to_send_to_ranks[p_id].push_back(varia);
            //====================Hybrid=======================
            for(int k=0;k<nvPerEl;k++)//This works for hexes.
            {
                v_id = ien->getVal(i,k);
                vertIDs_to_send_to_ranks[p_id].push_back(v_id);
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            for(int k=0;k<nfPerEl;k++)//This works for hexes.
            {
                f_id = ief->getVal(i,k);
                faceIDs_to_send_to_ranks[p_id].push_back(f_id);
            }
            //====================Hybrid=======================
            not_on_rank++;
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;
            nvPerEl = ie_Nv->getVal(i,0);
            nfPerEl = ie_Nf->getVal(i,0);
            for(int k=0;k<nvPerEl;k++)// looping over the vertices for element "i".
            {
                v_id = ien->getVal(i,k);
                
                //elem.push_back(v_id);
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_V_offsets,size,v_id);

                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vloc_tmp++;
                    }
                    lv_id++;
                }
            }
            for(int k=0;k<nfPerEl;k++)// looping over the vertices for element "i".
            {
                f_id = ief->getVal(i,k);
                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end() && f_id != -1) // add the required unique vertex for current rank.
                {
                    unique_faceIDs_on_rank_set.insert(f_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_F_offsets,size,f_id);

                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_face[r].push_back(f_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        faceIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        floc_tmp++;
                    }
                    lf_id++;
                }
            }
            
            loc_r_elem.push_back(el_id);
            loc_r_nv_elem.push_back(nvPerEl);
            loc_r_nf_elem.push_back(nfPerEl);
            loc_varia.push_back(varia);
            elem_set.insert(el_id);
            loc_r_elem_set.insert(el_id);
            elem_map[el_id] = on_rank;
            
            on_rank++;
        }
    }
    
    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);
    
    std::map<int,std::vector<int> >  part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> >  part_tot_recv_elNVs_map;
    std::map<int,std::vector<int> >  part_tot_recv_elNFs_map;
    
    std::map<int,std::vector<double> >  part_tot_recv_varias_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_f_map;
    std::map<int,std::vector<int> >::iterator it;
    
    
    
    
    int n_req_recv;
    int n_req_recv_v;
    int n_req_recv_f;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = elms_to_send_to_ranks.begin(); it != elms_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = vertIDs_to_send_to_ranks[it->first].size();
                int n_req_f         = faceIDs_to_send_to_ranks[it->first].size();
                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, comm);
                MPI_Send(&n_req_f, 1, MPI_INT, dest, dest*222, comm);

                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&varia_to_send_to_ranks[it->first][0], n_req, MPI_DOUBLE, dest, 20000+100+dest*2, comm);
                MPI_Send(&nvPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*33333+7777, comm);
                MPI_Send(&nfPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*44444+8888, comm);

                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222, comm, MPI_STATUS_IGNORE);

            std::vector<double> part_recv_varia(n_req_recv);
            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_el_nv(n_req_recv);
            std::vector<int>    part_recv_el_nf(n_req_recv);

            std::vector<int>    part_recv_vrt_id(n_req_recv_v);
            std::vector<int>    part_recv_face_id(n_req_recv_f);
            
            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_varia[0], n_req_recv, MPI_DOUBLE, q, 20000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_nv[0], n_req_recv, MPI_INT, q, rank*33333+7777, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_nf[0], n_req_recv, MPI_INT, q, rank*44444+8888, comm, MPI_STATUS_IGNORE);

            MPI_Recv(&part_recv_vrt_id[0],  n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            TotRecvElement_IDs_v_map[q] = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q] = part_recv_face_id;
            part_tot_recv_elIDs_map[q]  = part_recv_el_id;
            part_tot_recv_varias_map[q] = part_recv_varia;
            part_tot_recv_elNVs_map[q]  = part_recv_el_nv;
            part_tot_recv_elNFs_map[q]  = part_recv_el_nf;
        }
    }
    
    
    
    
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvElement_NVs;
    std::vector<int> TotRecvElement_NFs;
    std::vector<int> TotRecvVerts_IDs;
    std::vector<int> TotRecvFaces_IDs;
    std::vector<double> TotRecvElement_varia;
    
    std::map<int,std::vector<int> >::iterator totrecv;
    //unpack the element IDs and their corresponding variable values.
    int TotNelem_recv = 0;
    for(totrecv=part_tot_recv_elIDs_map.begin();totrecv!=part_tot_recv_elIDs_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElement_IDs.push_back(part_tot_recv_elIDs_map[totrecv->first][r]);
            TotRecvElement_NVs.push_back(part_tot_recv_elNVs_map[totrecv->first][r]);
            TotRecvElement_NFs.push_back(part_tot_recv_elNFs_map[totrecv->first][r]);
            TotRecvElement_varia.push_back(part_tot_recv_varias_map[totrecv->first][r]);
        }
        TotNelem_recv = TotNelem_recv + totrecv->second.size();
    }
    //unpack the vertex IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_v_map.begin();totrecv!=TotRecvElement_IDs_v_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvVerts_IDs.push_back(TotRecvElement_IDs_v_map[totrecv->first][r]);
        }
    }
    //unpack the face IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_f_map.begin();totrecv!=TotRecvElement_IDs_f_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvFaces_IDs.push_back(TotRecvElement_IDs_f_map[totrecv->first][r]);
        }
    }
    
    int Nel_extra = TotNelem_recv;
    int cnt_v = 0;
    int cnt_f = 0;
    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;
        int nvPerEl = TotRecvElement_NVs[i];
        int nfPerEl = TotRecvElement_NFs[i];
        for(int k=0;k<nvPerEl;k++)
        {
            int v_id_n = TotRecvVerts_IDs[cnt_v+k];
            //elem.push_back(v_id_n);
            r = FindRank(new_V_offsets,size,v_id_n);
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_V_offsets,size,v_id_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_vert[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    vloc_tmp++;
                }
                lv_id++;
            }
        }
        for(int k=0;k<nfPerEl;k++)// looping over the vertices for element "i".
        {
            int f_id_n = TotRecvVerts_IDs[cnt_v+k];
            if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_faceIDs_on_rank_set.insert(f_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_F_offsets,size,f_id_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_face[r].push_back(f_id_n); // add the vertex id that needs to be requested from rank r.
                }
                else
                {
                    faceIDs_on_rank.push_back(f_id_n);  // add the vertex to list that is already available on rank.
                    floc_tmp++;
                }
                lf_id++;
            }
        }
        
        cnt_v=cnt_v+nvPerEl;
        cnt_f=cnt_f+nfPerEl;
        
        //part_elem2verts.push_back(elem);
        elem_set.insert(TotRecvElement_IDs[i]);
        loc_r_elem_set.insert(TotRecvElement_IDs[i]);
        elem_map[el_id] = on_rank;
        on_rank++;
        
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
    int m = 0;
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);
    
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
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
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
    

    
    
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                double* vert_send = new double[nv_send*3];
                offset_xcn        = xcn_pstate->getOffset(rank);
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
                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
                delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
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

    int gvid=0;
    int lvid=0;


    for(m=0;m<vloc_tmp;m++)
    {
        gvid = vertIDs_on_rank[m];
       
        Vert* V = new Vert;
        
        V->x = xcn->getVal(gvid-xcn_o,0);
        V->y = xcn->getVal(gvid-xcn_o,1);
        V->z = xcn->getVal(gvid-xcn_o,2);
        
        LocalVerts.push_back(V);
        LocalVert2GlobalVert[lvid] = gvid;
        GlobalVert2LocalVert[gvid] = lvid;
        lvid++;
    }
    
    int o = 3*vloc_tmp;
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = rank2req_vert[it_f->first][u];
            
            Vert* V = new Vert;
            
            V->x = it_f->second[u*3+0];
            V->y = it_f->second[u*3+1];
            V->z = it_f->second[u*3+2];
            
            LocalVerts.push_back(V);
            
            LocalVert2GlobalVert[lvid]=gvid;
            GlobalVert2LocalVert[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }


    nLoc_Verts = LocalVerts.size();
    // ================================== Faces on Rank =========================================
    
    int lfid = 0;
    int gfid = 0;
    for(m=0;m<floc_tmp;m++)
    {
        gfid = faceIDs_on_rank[m];
    
        LocalFace2GlobalFace[lfid] = gfid;
        GlobalFace2LocalFace[gfid] = lfid;
        lfid++;
    }
    
    // ================================== Faces on Rank =========================================
    //NlocElem             = loc_r_elem.size()+Nel_extra+itel;
    //LocalElem2GlobalVert = new Array<int>(NlocElem,8);
    //LocalElem2LocalVert  = new Array<int>(NlocElem,8);
    //std::vector<double> U0vert;
    //U0Elem               = new Array<double>(NlocElem,1);
    //U0Vert               = new Array<double>(LocalVerts.size(),1);
    //ElemPart             = new Array<int>(NlocElem,1);

    int glob_v = 0;
    int loc_v  = 0;
    int glob_f = 0;
    int loc_f  = 0;
    double varia_v = 0.0;
    
    std::vector<int> tmp_globv;
    std::vector<int> tmp_locv;
    
    

    
    for(m=0;m<loc_r_elem.size();m++)
    {
        el_id   = loc_r_elem[m];
        nvPerEl = loc_r_nv_elem[m];
        nfPerEl = loc_r_nf_elem[m];
        varia_v = loc_varia[m];
        
        Loc_Elem.push_back(el_id);
        LocElem2Nv[el_id] = nvPerEl;
        LocElem2Nf[el_id] = nfPerEl;
        Loc_Elem_Nv.push_back(nvPerEl);
        Loc_Elem_Nf.push_back(nfPerEl);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(nvPerEl);
        LocAndAdj_Elem_Nf.push_back(nfPerEl);
        Loc_Elem_Varia.push_back(varia_v);
        LocAndAdj_Elem_Varia.push_back(varia_v);
        //U0Elem->setVal(m,0,rho_v);
        //ElemPart->setVal(m,0,el_id);
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;

        for(int p=0;p<nvPerEl;p++)
        {
            //GlobalFace2LocalFace
            //LocalVert2GlobalVert
            glob_v = ien->getVal(el_id-ien_o,p);
            if(glob_v>xcn_pstate->getNel())
            {
                std::cout << "Nel On error " << glob_v << std::endl;
            }
            
            loc_v  = GlobalVert2LocalVert[glob_v];
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            //LocalElem2GlobalVert->setVal(m,p,glob_v);
            //LocalElem2LocalVert->setVal(m,p,loc_v);
            //collect_var[loc_v].push_back(rho_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globVerts2globElem[glob_v].push_back(el_id);
            globElem2locVerts[el_id].push_back(loc_v);
        }
        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = ief->getVal(el_id-ien_o,p);
            loc_f  = GlobalFace2LocalFace[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
        }
        
        LocalElem2GlobalVert.push_back(tmp_globv);
        LocalElem2LocalVert.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }
    int cnv = 0;
    int cnf = 0;
    for(m=0;m<Nel_extra;m++)
    {
        el_id = TotRecvElement_IDs[m];
        nvPerEl = TotRecvElement_NVs[m];
        nfPerEl = TotRecvElement_NFs[m];
        LocElem2Nv[el_id] = nvPerEl;
        LocElem2Nf[el_id] = nfPerEl;
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        varia_v = TotRecvElement_varia[m];
        Loc_Elem.push_back(el_id);
        Loc_Elem_Nv.push_back(nvPerEl);
        Loc_Elem_Nf.push_back(nfPerEl);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(nvPerEl);
        LocAndAdj_Elem_Nf.push_back(nfPerEl);
        Loc_Elem_Varia.push_back(varia_v);
        LocAndAdj_Elem_Varia.push_back(varia_v);
        for(int p=0;p<nvPerEl;p++)
        {
            glob_v = TotRecvVerts_IDs[cnv+p];
            if(glob_v>xcn_pstate->getNel())
            {
                std::cout << "Nel Extra error " << glob_v << std::endl;
            }
            loc_v = GlobalVert2LocalVert[glob_v];
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            globVerts2globElem[glob_v].push_back(el_id);

            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            
        }
        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = TotRecvFaces_IDs[cnf+p];
            
            loc_f  = GlobalFace2LocalFace[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
            
        }
        cnv=cnv+nvPerEl;
        cnf=cnf+nfPerEl;
        LocalElem2GlobalVert.push_back(tmp_globv);
        LocalElem2LocalVert.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }
    
    nLoc_Elem = Loc_Elem.size();
    vloc = LocalVerts.size();
    floc = cnf;
    
    
    delete part_schedule;
    part_schedule_elem;
    
    TotRecvElement_IDs.clear();
    part_tot_recv_elIDs_map.clear();
    part_tot_recv_elNVs_map.clear();
    part_tot_recv_elNFs_map.clear();
    part_tot_recv_varias_map.clear();
    TotRecvElement_IDs_f_map.clear();
    TotRecvElement_IDs_v_map.clear();
    part_elem2verts.clear();
    vertIDs_to_send_to_ranks.clear();
    faceIDs_to_send_to_ranks.clear();
    varia_to_send_to_ranks.clear();
    rank2req_vert.clear();
    rank2req_face.clear();
    faceIDs_on_rank.clear();
    vertIDs_on_rank.clear();
    part_v.clear();
    delete[] new_V_offsets;
    delete[] new_F_offsets;
    TotRecvElement_IDs.clear();
    TotRecvElement_NVs.clear();
    TotRecvElement_NFs.clear();
    TotRecvVerts_IDs.clear();
    TotRecvFaces_IDs.clear();
    TotRecvElement_varia.clear();
    
    reqstd_ids_per_rank.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    tmp_globv.clear();
    tmp_locv.clear();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Partition::DetermineAdjacentElement2ProcMapUS3D(ParArray<int>* ien,
                                                     std::map<int,std::vector<int> > iee_vec,
                                                     ParArray<int>* part, ParArray<double>* xcn, Array<double>* U, MPI_Comm comm)
{
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int xcn_o = xcn_pstate->getOffset(rank);
    int ien_o = ien_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_face;
    int* new_V_offsets = new int[size];
    int* new_F_offsets = new int[size];
    int* new_E_offsets = new int[size];

    for(int i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
        new_F_offsets[i] = ife_pstate->getOffsets()[i]-1;
        new_E_offsets[i] = ien_pstate->getOffsets()[i]-1;
    }
        
    std::map<int,std::vector<int> > req_elem;
    
    int itel = 0;
    int Nel = part_global->getNrow();
    
    for(int i=0;i<Loc_Elem.size();i++)
    {
        int elId    = Loc_Elem[i];
        int nfPerEl = Loc_Elem_Nf[i];
        int k       = 0;
        
        for(int j=0;j<nfPerEl;j++)
        {
            int adjEl_id = iee_vec[elId][j];
            
            if((elem_set.find(adjEl_id)==elem_set.end()) && adjEl_id<Nel)
            {
                elem_set.insert(adjEl_id);
                p_id = part_global->getVal(adjEl_id,0);
                
                if(p_id != rank)
                {
                    adj_elements[p_id].push_back(adjEl_id);
                    req_elem[p_id].push_back(adjEl_id);
                    itel++;
                }
            }
        }
    }
    
    adj_schedule = DoScheduling(req_elem, comm);
    
    std::map<int,std::vector<int> >::iterator it;
    
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

                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req_adj_el, 1, MPI_INT, dest, 9876000+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req_adj_el, MPI_INT, dest, 9876000*2+dest*2, comm);
               i++;
            }
        }
        else if (adj_schedule->SendFromRank2Rank[q].find( rank ) != adj_schedule->SendFromRank2Rank[q].end())
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
    std::map<int,std::vector<int> > send_adj_NvertsPel;
    std::map<int,std::vector<int> > send_adj_NfacesPel;
    
    int TotNelem_adj_recv = 0;
    std::vector<int> TotAdj_El_IDs;
    int adj_id;

    int lelem = 0;

    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            adj_id = itv->second[j];
            int ladj_id = elem_map[adj_id];
            TotAdj_El_IDs.push_back(adj_id);

            int nvPerEl = LocElem2Nv[adj_id];
            int nfPerEl = LocElem2Nf[adj_id];
            send_adj_NvertsPel[dest].push_back(nvPerEl);
            send_adj_NfacesPel[dest].push_back(nfPerEl);
            for(int k=0;k<nvPerEl;k++)
            {
                v_id = globElem2globVerts[adj_id][k];
            
                send_adj_verts_IDs[dest].push_back(v_id);
            }
        
            for(int k=0;k<nfPerEl;k++)
            {
                f_id = globElem2globFaces[adj_id][k];
                send_adj_faces_IDs[dest].push_back(f_id);
            }
        }

        TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn   = 0;
    std::map<int,int  > recv_adj_back_Nverts;
    std::map<int,int* > recv_adj_back_verts_ids;
    std::map<int,int  > recv_adj_back_Nfaces;
    std::map<int,int* > recv_adj_back_faces_ids;
    //std::map<int,int > recv_adj_back_Nrhos;
    std::map<int,double* > recv_adj_back_rhos;
    std::map<int,std::vector<int>  > recv_adj_NvPel;
    std::map<int,std::vector<int>  > recv_adj_NfPel;
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
                
                int NnvPel = send_adj_NvertsPel[it->first].size();
                int NnfPel = send_adj_NvertsPel[it->first].size();
                MPI_Send(&NnvPel, 1, MPI_INT, dest, 98764444+5000*dest, comm);
                MPI_Send(&NnfPel, 1, MPI_INT, dest, 98764444-5000*dest, comm);
                
                MPI_Send(&send_adj_NvertsPel[it->first][0], NnvPel, MPI_INT, dest, 98364444+15000*dest, comm);
                MPI_Send(&send_adj_NfacesPel[it->first][0], NnfPel, MPI_INT, dest, 98364444-15000*dest, comm);
                
                int nf_adj_send = send_adj_faces_IDs[it->first].size();
                MPI_Send(&nf_adj_send, 1, MPI_INT, dest, 3333*9876+dest*8888,comm);
                MPI_Send(&send_adj_faces_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 2222*9876+dest*8888,comm);
                //int n_adj_rhos = send_adj_rhos[it->first].size();
                //MPI_Send(&n_adj_rhos, 1, MPI_INT, dest, 4444*9876+dest*8888,comm);
                //MPI_Send(&send_adj_rhos[it->first][0], n_adj_rhos, MPI_DOUBLE, dest, 5555*9876+dest*8888,comm);


            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 98760000+1000*rank, comm, MPI_STATUS_IGNORE);
            int* recv_adj_back_arr_ids = new int[n_adj_vert_recv_back];
            MPI_Recv(&recv_adj_back_arr_ids[0], n_adj_vert_recv_back, MPI_INT, q, 19999*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            
            int NnvPel_recv_back,NnfPel_recv_back;
            MPI_Recv(&NnvPel_recv_back, 1, MPI_INT, q, 98764444+5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NnfPel_recv_back, 1, MPI_INT, q, 98764444-5000*rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> Nnv_RB(NnvPel_recv_back);
            std::vector<int> Nnf_RB(NnfPel_recv_back);
            
            MPI_Recv(&Nnv_RB[0], NnvPel_recv_back, MPI_INT, q, 98364444+15000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&Nnf_RB[0], NnfPel_recv_back, MPI_INT, q, 98364444-15000*rank, comm, MPI_STATUS_IGNORE);
            
            int n_adj_face_recv_back;
            MPI_Recv(&n_adj_face_recv_back, 1, MPI_INT, q, 3333*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            int* recv_adj_back_arr_face_ids = new int[n_adj_face_recv_back];
            MPI_Recv(&recv_adj_back_arr_face_ids[0], n_adj_face_recv_back, MPI_INT, q, 2222*9876+rank*8888, comm,   MPI_STATUS_IGNORE);

            //int n_adj_rho_recv_back;
            //MPI_Recv(&n_adj_rho_recv_back, 1, MPI_INT, q, 4444*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            //double* recv_adj_back_arr_rho = new double[n_adj_rho_recv_back];
            //MPI_Recv(&recv_adj_back_arr_rho[0], n_adj_rho_recv_back, MPI_DOUBLE, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);


            recv_adj_back_Nverts[q]     = n_adj_vert_recv_back;
            recv_adj_back_verts_ids[q]  = recv_adj_back_arr_ids;

            recv_adj_back_Nfaces[q]     = n_adj_face_recv_back;
            recv_adj_back_faces_ids[q]  = recv_adj_back_arr_face_ids;

            recv_adj_NvPel[q] = Nnv_RB;
            recv_adj_NfPel[q] = Nnf_RB;
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
    std::vector<int> adj_elements_vec;
    std::vector<int> NvPEl_rb;
    std::vector<int> NfPEl_rb;
    std::map<int,std::vector<int> >::iterator itm_el;
    int offvvv = 0;
    for(itm_el=adj_elements.begin();itm_el!=adj_elements.end();itm_el++)
    {
        //TotNrho_adj_recv = TotNrho_adj_recv+itm->second;
        for(int i=0;i<itm_el->second.size();i++)
        {
            adj_elements_vec.push_back(adj_elements[itm_el->first][i]);
            int Nv = recv_adj_NvPel[itm_el->first][i];
            int Nf = recv_adj_NfPel[itm_el->first][i];
            NvPEl_rb.push_back(Nv);
            NfPEl_rb.push_back(Nf);
            offvvv=offvvv+Nv;

        }
    }

    //std::cout << "TotNelem_adj_recv " << TotNelem_adj_recv << " should be equal to " << TotNrho_adj_recv << std::endl;
    //std::cout << " Compare " <<  TotNelem_adj_recv << " " << itel << " " << adj_rhos.size() << std::endl;
    int cnt_v_adj = 0;
    int cnt_f_adj = 0;
    
    int offv = 0;
    int offf = 0;
    for(int i=0;i<itel;i++)
    {
        int Nv = NvPEl_rb[i];
        int Nf = NfPEl_rb[i];
        for(int k=0;k<Nv;k++)
        {
            //std::cout << offv+k << " " << adj_verts.size() << std::endl;
            int v_id_n = adj_verts[offv+k];
            
            r = FindRank(new_V_offsets,size,v_id_n);

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
                    vloc_tmp++;
                }
            }
        }

        for(int k=0;k<Nf;k++)
        {
            int f_id_n = adj_faces[offf+k];

            if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_faceIDs_on_rank_set.insert(f_id_n);
                faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                floc_tmp++;

                if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
                {
                    rank2req_face[r].push_back(f_id_n); // add vertex to rank2req_vert map.
                }
                else
                {
                    faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                    floc_tmp++;
                }
            }
        }
        offv = offv+Nv;
        offf = offf+Nf;
        //part_elem2verts.push_back(elem);
        //elem_set.insert(TotAdj_El_IDs[i]);
    }
    
    

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

       part_schedule = DoScheduling(rank2req_vert,comm);

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
                   MPI_Send(&n_req, 1, MPI_INT, dest, 6547+10*dest, comm);
                   //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                   MPI_Send(&it->second[0], n_req, MPI_INT, dest, 6547*2+dest*2, comm);

                   i++;
               }
           }
           else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
           {
               MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6547+10*rank, comm, MPI_STATUS_IGNORE);
               //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

               std::vector<int> recv_reqstd_ids(n_reqstd_ids);
               MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 6547*2+rank*2, comm, MPI_STATUS_IGNORE);
               reqstd_ids_per_rank[q] = recv_reqstd_ids;
           }
       }

       int offset_xcn = 0;
       int nloc_xcn = 0;
       std::map<int,int > recv_back_Nverts;
       std::map<int,double* > recv_back_verts;
       std::map<int,int* > recv_back_verts_ids;
       int n_recv_back;
    
    
    

    
    
       for(q=0;q<size;q++)
       {
           if(rank == q)
           {
               for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
               {
                   int nv_send = it->second.size();
                   double* vert_send = new double[nv_send*3];
                   offset_xcn        = xcn_pstate->getOffset(rank);
                   for(int u=0;u<it->second.size();u++)
                   {
                       vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                       vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                       vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);
                   }

                   int dest = it->first;
                   MPI_Send(&nv_send, 1, MPI_INT, dest, 6547+1000*dest, comm);
                   // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);

                   MPI_Send(&vert_send[0], nv_send*3, MPI_DOUBLE, dest, 6547+dest*8888, comm);
                   MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 8888*6547+dest*8888,comm);

                   delete[] vert_send;
               }
           }
           if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
            {
               MPI_Recv(&n_recv_back, 1, MPI_INT, q, 6547+1000*rank, comm, MPI_STATUS_IGNORE);

               double* recv_back_arr = new double[n_recv_back*3];
               int* recv_back_arr_ids = new int[n_recv_back];
               //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr[0], n_recv_back*3, MPI_DOUBLE, q, 6547+rank*8888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr_ids[0], n_recv_back, MPI_INT, q, 8888*6547+rank*8888, comm, MPI_STATUS_IGNORE);

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


       int gvid=0;
       int lvid=vloc;
       int m=0;
       for(m=0;m<vloc_tmp;m++)
       {
           gvid = vertIDs_on_rank[m];
           Vert* V = new Vert;

           V->x = xcn->getVal(gvid-xcn_o,0);
           V->y = xcn->getVal(gvid-xcn_o,1);
           V->z = xcn->getVal(gvid-xcn_o,2);

           LocalVerts.push_back(V);
           LocalVert2GlobalVert[lvid] = gvid;
           GlobalVert2LocalVert[gvid] = lvid;
           lvid++;
       }

       //int o = 3*vloc_tmp;
       m = 0;
       int u = 0;
       for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
       {
           int Nv = recv_back_Nverts[it_f->first];

           for(u=0;u<Nv;u++)
           {
               gvid = rank2req_vert[it_f->first][u];
               
               Vert* V = new Vert;

               V->x = it_f->second[u*3+0];
               V->y = it_f->second[u*3+1];
               V->z = it_f->second[u*3+2];

               LocalVerts.push_back(V);

               LocalVert2GlobalVert[lvid]=gvid;
               GlobalVert2LocalVert[gvid]=lvid;

               m++;
               lvid++;
           }
       }

       nLoc_Verts = LocalVerts.size();
    //std::cout << " " << rank << " LocalVerts.size() after " << LocalVerts.size() << std::endl;

    //NlocElem = eloc;
    //std::cout << "eloc " << eloc << std::endl;
    // ================================== Faces on Rank =========================================

    int lfid = floc;
    int gfid = 0;
    for(int m=0;m<floc_tmp;m++)
    {
        gfid = faceIDs_on_rank[m];

        LocalFace2GlobalFace[lfid] = gfid;
        GlobalFace2LocalFace[gfid] = lfid;
        lfid++;
    }

    int cnv = 0;
    int cnf = 0;
    int idsave = 0;
    //double rho_v;
    int loc_v;
    int glob_f;
    int loc_f;
    int glob_v;
    //std::cout << adj_verts.size() << " " << Nel_extra2 <<std::endl;
    std::vector<int> tmp_globv;
    std::vector<int> tmp_locv;
    std::vector<int> tmp_globf;
    std::vector<int> tmp_locf;
    

    
    
    double varia_v = 0.0;
    for(int m=0;m<itel;m++)
    {
        el_id = adj_elements_vec[m];
        int Nv = NvPEl_rb[m];
        int Nf = NfPEl_rb[m];
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        //rho_v = adj_rhos[m];

        //U0Elem.push_back(rho_v);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(Nv);
        LocAndAdj_Elem_Nf.push_back(Nf);
        //LocAndAdj_Elem_Varia.push_back(varia_v);
//      U0Elem->setVal(m+o,0,rho_v);
//      ElemPart->setVal(m+o,0,el_id);

        for(int p=0;p<Nv;p++)
        {
            glob_v = adj_verts[cnv];
            loc_v  = GlobalVert2LocalVert[glob_v];
            //LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            //LocalElem2LocalVert->setVal(m+o,p,loc_v);
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globVerts2globElem[glob_v].push_back(el_id);

            globElem2locVerts[el_id].push_back(loc_v);
            //collect_var[loc_v].push_back(rho_v);
            cnv++;
            
        }
        for(int p=0;p<Nf;p++)
        {
            glob_f = adj_faces[cnf];
            loc_f  = GlobalFace2LocalFace[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
            cnf++;
        }
        LocalElem2GlobalVert.push_back(tmp_globv);
        LocalElem2LocalVert.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }

//    nLoc_Elem = U0Elem.size();
//    std::map<int,std::vector<double> >::iterator it_rhos;
//    double sum = 0;
//    int c = 0;
//
//    for(it_rhos=collect_var.begin();it_rhos!=collect_var.end();it_rhos++)
//    {
//        sum = 0;
//        for(int q = 0;q<it_rhos->second.size();q++)
//        {
//            sum = sum + it_rhos->second[q];
//        }
//        U0Vert->setVal(c,0,sum/it_rhos->second.size());
//        c++;
//    }

    delete[] new_V_offsets;
    delete[] new_F_offsets;
    delete[] new_E_offsets;

    reqstd_ids_per_rank.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    
    NvPEl_rb.clear();
    NfPEl_rb.clear();
    
    rank2req_vert.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    tmp_globv.clear();
    tmp_locv.clear();
    tmp_globf.clear();
    tmp_locf.clear();

}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void Partition::AddStateForAdjacentElements(std::map<int,double> U, MPI_Comm comm)
{
    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //std::cout << xcn->getOffset(rank) << " " << xcn_pstate->getOffset(rank) << std::endl;
    
    //int ien_o = ien_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Elems;
    int* new_offsets = new int[size];
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ien_pstate->getOffsets()[i]-1;
    }
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    for(int i=0;i<LocAndAdj_Elem.size();i++)
    {
        int el_req = LocAndAdj_Elem[i];
        
        r = part_global->getVal(el_req,0);

        if(r != rank)
        {
            rank2req_Elems[r].push_back(el_req);
        }
    }

    ScheduleObj* iee_schedule = DoScheduling(rank2req_Elems,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Elems.begin(); it != rank2req_Elems.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+20*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*40, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+20*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*40, comm, MPI_STATUS_IGNORE);
            
            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    rank2req_Elems.clear();
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IEE_Elem_IDs;
    std::vector<int> TotIEE_El_IDs;


    int TotNelem_IEE_recv   = 0;
    int eIEE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Niee;
    std::map<int,int* > recv_back_el_ids;
    std::map<int,double* > recv_back_iee;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_E_IDs_per_rank.begin(); it != reqstd_E_IDs_per_rank.end(); it++)
            {
                int ne_send             = it->second.size();
                double* iee_send        = new double[ne_send];
                //int offset_iee          = ien_pstate->getOffset(rank);
                
                for(int u=0;u<it->second.size();u++)
                {
                    iee_send[u]=U[it->second[u]];
                }

                int dest = it->first;
                MPI_Send(&ne_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], ne_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&iee_send[0], ne_send, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                delete[] iee_send;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_iee_arr   = new double[n_recv_back];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_iee_arr[0], n_recv_back, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Niee[q]     = n_recv_back;
            recv_back_el_ids[q]   = recv_back_ids_arr;
            recv_back_iee[q]      = recv_back_iee_arr;

         }
    }
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Niee.begin();iter!=recv_back_Niee.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            el_id = recv_back_el_ids[iter->first][s];
            U[el_id] = recv_back_iee[iter->first][s];
        }
    }
    
//    recv_back_Niee.clear();
//    recv_back_iee.clear();
//    recv_back_el_ids.clear();
    delete[] new_offsets;
}



void Partition::AddStateVecForAdjacentElements(std::map<int,Array<double>* > &U, int nvar, MPI_Comm comm)
{
    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //std::cout << xcn->getOffset(rank) << " " << xcn_pstate->getOffset(rank) << std::endl;
    
    //int ien_o = ien_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Elems;
    int* new_offsets = new int[size];
    std::map<int,Array<double>* > U_loc;
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ien_pstate->getOffsets()[i]-1;
    }
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    for(int i=0;i<LocAndAdj_Elem.size();i++)
    {
        int el_req = LocAndAdj_Elem[i];
        
        r = part_global->getVal(el_req,0);

        if(r != rank)
        {
            rank2req_Elems[r].push_back(el_req);
        }
//        else
//        {
//            Array<double>* StateVec = new Array<double>(nvar,1);
//            for(int q=0;q<nvar;q++)
//            {
//                StateVec->setVal(q,0,U[el_req]->getVal(q,0));
//            }
//            U_loc[el_req] = StateVec;
//        }
    }
    
    ScheduleObj* iee_schedule = DoScheduling(rank2req_Elems,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Elems.begin(); it != rank2req_Elems.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+20*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*40, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+20*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*40, comm, MPI_STATUS_IGNORE);
            
            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IEE_Elem_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IEE_recv   = 0;
    int eIEE_id             = 0;
    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Niee;
    std::map<int,int* > recv_back_el_ids;
    std::map<int,double* > recv_back_iee;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_E_IDs_per_rank.begin(); it != reqstd_E_IDs_per_rank.end(); it++)
            {
                int ne_send             = it->second.size();
                double* iee_send        = new double[ne_send*nvar];
                //int offset_iee          = ien_pstate->getOffset(rank);
                
                for(int u=0;u<ne_send;u++)
                {
                    for(int s=0;s<nvar;s++)
                    {
                        iee_send[u*nvar+s]=U[it->second[u]]->getVal(s,0);
                    }
                }

                int dest = it->first;
                MPI_Send(&ne_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], ne_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&iee_send[0], ne_send*nvar, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                delete[] iee_send;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            int*    recv_back_ids_arr   = new int[n_recv_back];
            double* recv_back_iee_arr   = new double[n_recv_back*nvar];

            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_iee_arr[0], n_recv_back*nvar, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Niee[q]     = n_recv_back;
            recv_back_el_ids[q]   = recv_back_ids_arr;
            recv_back_iee[q]      = recv_back_iee_arr;

         }
    }
    
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Niee.begin();iter!=recv_back_Niee.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            el_id = recv_back_el_ids[iter->first][s];
            Array<double>* StateVec = new Array<double>(nvar,1);
            for(int p=0;p<nvar;p++)
            {
                StateVec->setVal(p,0,recv_back_iee[iter->first][s*nvar+p]);
            }
            U[el_id] = StateVec;
        }
    }
    
    delete[] new_offsets;
}


void Partition::AddStateVecForAdjacentVertices(std::map<int,Array<double>* > &Uv, int nvar, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> > send_adj_verts_IDs;
    std::map<int,std::vector<double> > send_adj_verts_Us;
    
    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)// private
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            int adj_id  = itv->second[j];
            //int ladj_id = elem_map[adj_id]; // private

            int nvPerEl = LocElem2Nv[adj_id];
            
            for(int k=0;k<nvPerEl;k++)
            {
                int v_id = globElem2globVerts[adj_id][k];
                send_adj_verts_IDs[dest].push_back(v_id);
                for(int l=0;l<nvar;l++)
                {
                    send_adj_verts_Us[dest].push_back(Uv[v_id]->getVal(l,0));
                }
            }
        }
        //TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn   = 0;
    int n_adj_vert_recv_back;
    std::map<int,std::vector<int> > recv_adj_back_verts_ids;
    std::map<int,std::vector<double> > recv_adj_back_verts_Us;
    
    // This sends the right vertices of the requested elements to correct processor.
    
    
    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = send_adj_verts_IDs.begin(); it != send_adj_verts_IDs.end(); it++)
            {
                int nv_adj_send       = it->second.size();
                int dest = it->first;
                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 12340000+1000*dest, comm);
                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 19999*1234+dest*8888,comm);
                MPI_Send(&send_adj_verts_Us[it->first][0], it->second.size()*nvar, MPI_DOUBLE, dest, 1234*1234+dest*8888,comm);
            
                


            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 12340000+1000*rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<int>recv_adj_back_vec_ids(n_adj_vert_recv_back);
            MPI_Recv(&recv_adj_back_vec_ids[0], n_adj_vert_recv_back, MPI_INT, q, 19999*1234+rank*8888, comm, MPI_STATUS_IGNORE);
            
            std::vector<double>recv_adj_back_vec_Us(n_adj_vert_recv_back*nvar);
            MPI_Recv(&recv_adj_back_vec_Us[0], n_adj_vert_recv_back*nvar, MPI_DOUBLE, q, 1234*1234+rank*8888, comm, MPI_STATUS_IGNORE);
                        
            recv_adj_back_verts_ids[q] = recv_adj_back_vec_ids;
            recv_adj_back_verts_Us[q]  = recv_adj_back_vec_Us;

        }
    }
    
    
    int TotNvert_adj_recv = 0;
    
    std::map<int,std::vector<int> >::iterator itm;
    for(itm=recv_adj_back_verts_ids.begin();itm!=recv_adj_back_verts_ids.end();itm++)
    {
        TotNvert_adj_recv = TotNvert_adj_recv+itm->second.size();
        for(int i=0;i<itm->second.size();i++)
        {
            if(Uv.find(recv_adj_back_verts_ids[itm->first][i])==Uv.end())
            {
                Array<double>* StateVec = new Array<double>(nvar,1);
                
                for(int u=0;u<nvar;u++)
                {
                    StateVec->setVal(u,0,recv_adj_back_verts_Us[itm->first][i*nvar+u]);
                }
                
                Uv[recv_adj_back_verts_ids[itm->first][i]] = StateVec;
            }
        }
    }
    
    send_adj_verts_Us.clear();
    recv_adj_back_verts_ids.clear();
    recv_adj_back_verts_Us.clear();
    send_adj_verts_IDs.clear();

}



void Partition::AddAdjacentVertexDataUS3D(std::map<int,double> &Uv, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,double> Uv_adj = Uv;
    std::map<int,std::vector<double> > send_adj_aux;
    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> > send_adj_verts_IDs;
    std::map<int,std::vector<double> > send_adj_verts_Us;
    

    
    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)// private
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            int adj_id  = itv->second[j];
            //int ladj_id = elem_map[adj_id]; // private

            int nvPerEl = LocElem2Nv[adj_id];
            //int nfPerEl = LocElem2Nf[adj_id];
            
            for(int k=0;k<nvPerEl;k++)
            {
                int v_id = globElem2globVerts[adj_id][k];
                send_adj_verts_IDs[dest].push_back(v_id);
                send_adj_verts_Us[dest].push_back(Uv[v_id]);
            }
        }
    }
    
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn   = 0;
    int n_adj_vert_recv_back;
    std::map<int,std::vector<int> > recv_adj_back_verts_ids;
    std::map<int,std::vector<double> > recv_adj_back_verts_Us;
    
    // This sends the right vertices of the requested elements to correct processor.
    
    
    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = send_adj_verts_IDs.begin(); it != send_adj_verts_IDs.end(); it++)
            {
                int nv_adj_send       = it->second.size();
                int dest = it->first;
                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 12340000+1000*dest, comm);
                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 19999*1234+dest*8888,comm);
                MPI_Send(&send_adj_verts_Us[it->first][0], it->second.size(), MPI_DOUBLE, dest, 1234*1234+dest*8888,comm);
            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 12340000+1000*rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<int>recv_adj_back_vec_ids(n_adj_vert_recv_back);
            MPI_Recv(&recv_adj_back_vec_ids[0], n_adj_vert_recv_back, MPI_INT, q, 19999*1234+rank*8888, comm, MPI_STATUS_IGNORE);
            
            std::vector<double>recv_adj_back_vec_Us(n_adj_vert_recv_back);
            MPI_Recv(&recv_adj_back_vec_Us[0], n_adj_vert_recv_back, MPI_DOUBLE, q, 1234*1234+rank*8888, comm, MPI_STATUS_IGNORE);
                        
            recv_adj_back_verts_ids[q] = recv_adj_back_vec_ids;
            recv_adj_back_verts_Us[q]  = recv_adj_back_vec_Us;

        }
    }
    
    int TotNvert_adj_recv = 0;
    
    std::map<int,std::vector<int> >::iterator itm;
    for(itm=recv_adj_back_verts_ids.begin();itm!=recv_adj_back_verts_ids.end();itm++)
    {
        TotNvert_adj_recv = TotNvert_adj_recv+itm->second.size();
        for(int i=0;i<itm->second.size();i++)
        {
            if(Uv.find(recv_adj_back_verts_ids[itm->first][i])==Uv.end())
            {
                Uv[recv_adj_back_verts_ids[itm->first][i]] = recv_adj_back_verts_Us[itm->first][i];
            }
        }
    }
    
    
    
    send_adj_aux.clear();

    
    send_adj_verts_Us.clear();
    recv_adj_back_verts_ids.clear();
    recv_adj_back_verts_Us.clear();
    send_adj_verts_IDs.clear();
    
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



std::vector<double> Partition::PartitionAuxilaryData(Array<double>* U, MPI_Comm comm)
{
    
    std::vector<double> UauxElem;
    int i;
    // First send the aux data based on the partition array/
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<double> > aux_to_send_to_ranks;
    std::map<int,double> elid_2_auxVal;
    std::vector<double> aux_on_rank;
    

    double aux = 0.0;
    int p_id;
    int el_id;
    for(i=0;i<part->getNrow();i++)
    {
        p_id  = part->getVal(i,0);
        el_id = ien_pstate->getOffset(rank)+i;
        aux   = U->getVal(i,0);
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            aux_to_send_to_ranks[p_id].push_back(aux);
        }
        
        UauxElem.push_back(aux);
        elid_2_auxVal[el_id] = aux;


    }
    
    ScheduleObj* aux_schedule = DoScheduling(elms_to_send_to_ranks, comm);

    std::map<int,std::vector<double> > recv_FromRanks_aux;
    std::map<int,std::vector<double> >::iterator it;
    int n_req_recv;
    int n_req_recv_v;
    
    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = aux_to_send_to_ranks.begin(); it != aux_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                                
                MPI_Send(&n_req, 1, MPI_INT, dest, dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_DOUBLE, dest, 20000+100+dest*2, comm);
                
            }
        }
        else if (aux_schedule->SendFromRank2Rank[q].find( rank ) != aux_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<double> recv_aux(n_req_recv);
            MPI_Recv(&recv_aux[0], n_req_recv, MPI_DOUBLE, q, 20000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            recv_FromRanks_aux[q] = recv_aux;
        }
    }
    
    std::map<int,std::vector<double> >::iterator it2;
    for(it2 = recv_FromRanks_aux.begin();it2!=recv_FromRanks_aux.end();it2++)
    {
        for(int s=0;s<it2->second.size();s++)
        {
            el_id = part_tot_recv_elIDs[it2->first][s];
            UauxElem.push_back(it2->second[s]);
            elid_2_auxVal[el_id] = it2->second[s];
        }
    }
    
    // Once the aux data is send and received by all processors based on the partition array.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // reqstd_adj_ids_per_rank is data structure that is part of the Partition Object;
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::map<int,std::vector<double> > send_adj_aux;
    int adj_id;
    std::map<int,std::vector<int> >::iterator itv;
    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            adj_id = itv->second[j];
            send_adj_aux[dest].push_back(elid_2_auxVal[adj_id]);
        }
    }
    
    // This sends the right vertices of the requested elements to correct processor.
    std::map<int,std::vector<double> > recv_adj_back_aux;
    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = send_adj_aux.begin(); it != send_adj_aux.end(); it++)
            {
                int nv_adj_send       = it->second.size();
                int dest = it->first;
                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 987600+1000*dest, comm);
                MPI_Send(&it->second[0], it->second.size(), MPI_DOUBLE, dest, 17777*9876+dest*8888,comm);
            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            int n_adj_aux_recv_back;
            MPI_Recv(&n_adj_aux_recv_back, 1, MPI_INT, q, 987600+1000*rank, comm, MPI_STATUS_IGNORE);
            std::vector<double> recv_adj_back_vec_aux(n_adj_aux_recv_back);
            MPI_Recv(&recv_adj_back_vec_aux[0], n_adj_aux_recv_back, MPI_DOUBLE, q, 17777*9876+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_adj_back_aux[q]  = recv_adj_back_vec_aux;

        }
    }
    
    for(it2 = recv_adj_back_aux.begin();it2!=recv_adj_back_aux.end();it2++)
    {
        for(int s=0;s<it2->second.size();s++)
        {
            el_id = adj_elements[it2->first][s];
            UauxElem.push_back(it2->second[s]);
            elid_2_auxVal[el_id] = aux;
        }
    }
    
    /*
    recv_FromRanks_aux.clear();
    recv_adj_back_aux.clear();
    elid_2_auxVal.clear();
    aux_to_send_to_ranks.clear();
    aux_on_rank.clear();
    */
    
    return UauxElem;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



i_part_map* Partition::getElement2EntityPerPartition(ParArray<int>* iee, std::vector<int> Loc_Elem_Ne, MPI_Comm comm)
{
    
    i_part_map* iee_p_map = new i_part_map;
    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    //std::cout << xcn->getOffset(rank) << " " << xcn_pstate->getOffset(rank) << std::endl;

    int ien_o = ien_pstate->getOffset(rank);
    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Elems;
    std::map<int,std::vector<int> > rank2req_Elems_Ne;
    int* new_offsets = new int[size];
    std::map<int,std::vector<int> > iee_loc;
    std::map<int,std::vector<int> > iee_loc_inv;
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ien_pstate->getOffsets()[i]-1;
    }
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    for(int i=0;i<Loc_Elem.size();i++)
    {
        int el_req          = Loc_Elem[i];
        int nEntityPelement = Loc_Elem_Ne[i];
        r                   = FindRank(new_offsets,size,el_req);
        
        if(r != rank)
        {
            rank2req_Elems[r].push_back(el_req);
            rank2req_Elems_Ne[r].push_back(nEntityPelement);
        }
        else
        {
            for(int j=0;j<nEntityPelement;j++)
            {
                iee_loc[el_req].push_back(iee->getVal(el_req-ien_o,j));
                iee_loc_inv[iee->getVal(el_req-ien_o,j)].push_back(el_req);
            }
        }
    }

    ScheduleObj* iee_schedule = DoScheduling(rank2req_Elems,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;
    std::map<int,std::vector<int> >  reqstd_E_NePID_per_rank;
    
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Elems.begin(); it != rank2req_Elems.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);
                MPI_Send(&rank2req_Elems_Ne[it->first][0], n_req, MPI_INT, dest, 6611+dest*2, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_NePid(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_reqstd_NePid[0], n_reqstd_ids, MPI_INT, q, 6611+rank*2, comm, MPI_STATUS_IGNORE);

            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
            reqstd_E_NePID_per_rank[q] = recv_reqstd_NePid;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IEE_recv   = 0;
    int eIEE_id             = 0;
    
    int offset_xcn          = 0;
    int nloc_xcn            = 0;
    
    std::map<int,std::vector<int> > recv_back_el_ids;
    std::map<int,std::vector<int> > recv_back_el_NePids;
    std::map<int,std::vector<int> > recv_back_iee;
    
    int n_recv_back;
    int nNe_recv_back;
    int offs = 0;
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_E_IDs_per_rank.begin(); it != reqstd_E_IDs_per_rank.end(); it++)
            {
                int dest = it->first;

                int ne_send       = it->second.size();
                int offset_iee    = ien_pstate->getOffset(rank);
                int nePid_t       = 0;
                for(int u=0;u<ne_send;u++)
                {
                    int ne_p_id = reqstd_E_NePID_per_rank[it->first][u];
                    nePid_t = nePid_t+ne_p_id;
                }
                
                int* iee_send  = new int[nePid_t];
                offs = 0;
                for(int u=0;u<ne_send;u++)
                {
                    int ncol = reqstd_E_NePID_per_rank[it->first][u];
                    for(int h=0;h<ncol;h++)
                    {
                        iee_send[offs+h] = iee->getVal(it->second[u]-offset_iee,h);
                    }
                    offs=offs+ncol;
                }
                
                MPI_Send(&ne_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                MPI_Send(&nePid_t, 1, MPI_INT, dest, 9876*2222+1000*dest, comm);
                
                MPI_Send(&it->second[0], ne_send, MPI_INT, dest, 9876*7777+dest*888, comm);
                MPI_Send(&reqstd_E_NePID_per_rank[it->first][0], ne_send, MPI_INT, dest, 9876*2222+dest*1000, comm);
                MPI_Send(&iee_send[0], nePid_t, MPI_INT, dest, 9876*6666+dest*8888, comm);

                //delete[] iee_send;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nNe_recv_back, 1, MPI_INT, q, 9876*2222+1000*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_back_ids_arr(n_recv_back);
            std::vector<int> recv_back_NePids_arr(n_recv_back);
             
            std::vector<int> recv_back_iee_arr(nNe_recv_back);


            MPI_Recv(&recv_back_ids_arr[0],     n_recv_back,      MPI_INT, q, 9876*7777+rank*888,  comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_NePids_arr[0],  n_recv_back,      MPI_INT, q, 9876*2222+rank*1000, comm, MPI_STATUS_IGNORE);
             
            MPI_Recv(&recv_back_iee_arr[0],     nNe_recv_back,    MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_el_ids[q]     = recv_back_ids_arr;
            recv_back_el_NePids[q]  = recv_back_NePids_arr;
            recv_back_iee[q]        = recv_back_iee_arr;
         }
    }

    std::map<int,std::vector<int> >::iterator iter;
    int ntotal=0;
    ee.clear();
    
    for(iter=recv_back_el_ids.begin();iter!=recv_back_el_ids.end();iter++)
    {
        int L = iter->second.size();
        int offs = 0;
        for(int s=0;s<L;s++)
        {
            el_id = iter->second[s];
            int NePid = recv_back_el_NePids[iter->first][s];
            for(int r=0;r<NePid;r++)
            {
                iee_loc[el_id].push_back(recv_back_iee[iter->first][offs+r]);
                iee_loc_inv[recv_back_iee[iter->first][offs+r]].push_back(el_id);
    
            }
            offs = offs+NePid;
        }

        
        ntotal=ntotal+L;
    }
    delete[] new_offsets;
    
    iee_p_map->i_map = iee_loc;
    iee_p_map->i_inv_map = iee_loc_inv;
    
    return iee_p_map;
}




i_part_map* Partition::getFace2EntityPerPartition(ParArray<int>* ife,MPI_Comm comm)
{
    
    i_part_map* ife_p_map = new i_part_map;
    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int ife_o = ife_pstate->getOffset(rank);
    int face_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Faces;
    int* new_offsets = new int[size];
    std::map<int,std::vector<int> > ife_loc;
    std::map<int,std::vector<int> > ife_loc_inv;
    
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    int ncol = ife->getNcol();
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;
    
    for(itefmap=ief_part_map->i_map.begin();itefmap!=ief_part_map->i_map.end();itefmap++)
    {
        for(int q=0;q<itefmap->second.size();q++)
        {
            int face_req = itefmap->second[q];
            
            r = FindRank(new_offsets,size,face_req);
            
            if(r != rank)
            {
                rank2req_Faces[r].push_back(face_req);
            }
            else
            {
                for(int j=0;j<ncol;j++)
                {
                    ife_loc[face_req].push_back(ife->getVal(face_req-ife_o,j));
                    ife_loc_inv[ife->getVal(face_req-ife_o,j)].push_back(face_req);
                }
            }
        }
    }
    
    
    
    
    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nife;
    std::map<int,int* > recv_back_face_ids;
    std::map<int,double* > recv_back_ife;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send             = it->second.size();
                double* ife_send        = new double[nf_send*ncol];
                int offset_ife          = ife_pstate->getOffset(rank);
                
                for(int u=0;u<it->second.size();u++)
                {
                    for(int h=0;h<ncol;h++)
                    {
                        ife_send[u*ncol+h] = ife->getVal(it->second[u]-offset_ife,h);
                    }
                }

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&ife_send[0], nf_send*ncol, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                delete[] ife_send;
            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_ife_arr   = new double[n_recv_back*ncol];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_ife_arr[0], n_recv_back*ncol, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nife[q]       = n_recv_back;
            recv_back_face_ids[q]   = recv_back_ids_arr;
            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Nife.begin();iter!=recv_back_Nife.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            face_id = recv_back_face_ids[iter->first][s];
            for(int r=0;r<ncol;r++)
            {
                ife_loc[face_id].push_back(recv_back_ife[iter->first][s*ncol+r]);
                ife_loc_inv[recv_back_ife[iter->first][s*ncol+r]].push_back(face_id);
                
            }
        }
        ntotal=ntotal+L;
    }

    delete[] new_offsets;
    
    ife_p_map->i_map = ife_loc;
    ife_p_map->i_inv_map = ife_loc_inv;
    
    return ife_p_map;
}




i_part_map* Partition::getFace2NodePerPartition(ParArray<int>* ifn, MPI_Comm comm)
{
    
    i_part_map* ife_p_map = new i_part_map;

    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int ife_o = ife_pstate->getOffset(rank);
    int face_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Faces;
    std::map<int,std::vector<int> > rank2req_FacesNv;
    int* new_offsets = new int[size];
    std::map<int,std::vector<int> > ife_loc;
    std::map<int,std::vector<int> > ife_loc_inv;

    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }

    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;

    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;

    for(itefmap=ief_part_map->i_map.begin();itefmap!=ief_part_map->i_map.end();itefmap++)
    {
        
        for(int q=0;q<itefmap->second.size();q++)
        {
            int face_req = itefmap->second[q];
            int Nv = if_Nv_part_map->i_map[face_req][0];
            r = FindRank(new_offsets,size,face_req);

            if(r != rank)
            {
                rank2req_Faces[r].push_back(face_req);
                rank2req_FacesNv[r].push_back(Nv);
            }
            else
            {
                for(int j=0;j<Nv;j++)
                {
                    ife_loc[face_req].push_back(ifn->getVal(face_req-ife_o,j));
                    ife_loc_inv[ifn->getVal(face_req-ife_o,j)].push_back(face_req);
                }
            }
        }
    }
    
    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;
    std::map<int,std::vector<int> >  reqstd_F_Nvs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;


                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);
                MPI_Send(&rank2req_FacesNv[it->first][0], n_req, MPI_INT, dest, 3876*2*7654+dest*2, comm);
                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_Nvs(n_reqstd_ids);

            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_reqstd_Nvs[0], n_reqstd_ids, MPI_INT, q, 3876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);

            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
            reqstd_F_Nvs_per_rank[q] = recv_reqstd_Nvs;
        }
    }

    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,std::vector<int> > recv_back_fNvs;
    std::map<int,std::vector<int> > recv_back_fids;
    std::map<int,std::vector<int> > recv_back_ife;
    int n_recv_back;
    int nvsPface;
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send             = it->second.size();
                int offset_ife          = ife_pstate->getOffset(rank);
                int nvPf_t       = 0;
                for(int u=0;u<nf_send;u++)
                {
                    nvPf_t = nvPf_t + reqstd_F_Nvs_per_rank[it->first][u];
                }
                
                int* ifn_send        = new int[nvPf_t];
                int offs = 0;
                for(int u=0;u<it->second.size();u++)
                {
                    int ncol = reqstd_F_Nvs_per_rank[it->first][u];
                    for(int h=0;h<ncol;h++)
                    {
                        ifn_send[offs+h] = ifn->getVal(it->second[u]-offset_ife,h);
                    }
                    offs=offs+ncol;
                }

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                MPI_Send(&nvPf_t, 1, MPI_INT, dest, 9876*33+1000*dest, comm);

                MPI_Send(&it->second[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);
                
                MPI_Send(&reqstd_F_Nvs_per_rank[it->first][0], nf_send, MPI_INT, dest, 9876*1111+dest*888, comm);

                MPI_Send(&ifn_send[0], nvPf_t, MPI_INT, dest, 9876*6666+dest*8888, comm);

                delete[] ifn_send;
            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nvsPface, 1, MPI_INT, q, 9876*33+1000*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_back_ids_vec(n_recv_back);
            std::vector<int> recv_back_Nvs_vec(n_recv_back);
            std::vector<int> recv_back_ife_vec(nvsPface);

            MPI_Recv(&recv_back_ids_vec[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
             
            MPI_Recv(&recv_back_Nvs_vec[0], n_recv_back, MPI_INT, q, 9876*1111+rank*888, comm, MPI_STATUS_IGNORE);
            
            MPI_Recv(&recv_back_ife_vec[0], nvsPface, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_fNvs[q]       = recv_back_Nvs_vec;
            recv_back_fids[q]       = recv_back_ids_vec;
            recv_back_ife[q]        = recv_back_ife_vec;

         }
    }
//
    
    std::map<int,std::vector<int> >::iterator iter;
    ee.clear();
    for(iter=recv_back_fids.begin();iter!=recv_back_fids.end();iter++)
    {
        int L = iter->second.size();
        int offss = 0;
        for(int s=0;s<L;s++)
        {
            face_id  = recv_back_fids[iter->first][s];
            int ncol = recv_back_fNvs[iter->first][s];
            
            for(int r=0;r<ncol;r++)
            {
                ife_loc[face_id].push_back(recv_back_ife[iter->first][offss+r]);
                ife_loc_inv[recv_back_ife[iter->first][offss+r]].push_back(face_id);
            }
            offss = offss+ncol;
        }
    }

    delete[] new_offsets;

    ife_p_map->i_map = ife_loc;
    ife_p_map->i_inv_map = ife_loc_inv;
    
    return ife_p_map;
}

void Partition::CreatePartitionDomain()
{
    pDom = new Domain;
    
    std::map<int,int> gv2lpv;
    std::map<int,int> lv2gpv;
    std::map<int,int> gv2lpartv;
    std::map<int,int> lpartv2gv;
    std::set<int> gv_set;
    Array<int>* locelem2locnode= new Array<int>(ien_part_map->i_map.size(),8);

    std::map<int,std::vector<int> >::iterator itm;
    std::map<int,std::vector<int> > vert2elem;
    
    int lcv  = 0;
    int elid = 0;
    
    std::vector<int> loc_part_verts;
    std::vector<int> glob_part_verts;
    
    std::map<int,std::vector<int> > Elements;
    std::map<int,std::vector<int> > Hexes;
    std::map<int,std::vector<int> > Prisms;
    std::map<int,std::vector<int> > Tetras;
    
    std::map<int,std::vector<int> > GElements;
    std::map<int,std::vector<int> > GHexes;
    std::map<int,std::vector<int> > GPrisms;
    std::map<int,std::vector<int> > GTetras;
    
    for(itm  = ien_part_map->i_map.begin();
        itm != ien_part_map->i_map.end();
        itm++)
    {
        int glob_id  = itm->first;
        std::vector<int>El(itm->second.size());
        std::vector<int>Elg(itm->second.size());
        for(int q=0;q<itm->second.size();q++)
        {
            int gv = itm->second[q];
            int lv = GlobalVert2LocalVert[gv];
            
            if(gv_set.find(gv)==gv_set.end())
            {
                gv_set.insert(gv);
                loc_part_verts.push_back(lv);
                vert2elem[gv].push_back(glob_id);
                gv2lpv[gv]=lcv;
                lv2gpv[lcv]=gv;
                gv2lpartv[gv]=lv;
                lpartv2gv[lv]=gv;
                locelem2locnode->setVal(elid,q,lcv);
                El[q] = lcv;
                Elg[q] = gv;
                lcv=lcv+1;
            }
            else
            {
                int lcv_u = gv2lpv[gv];
                locelem2locnode->setVal(elid,q,lcv_u);
                vert2elem[gv].push_back(glob_id);
                El[q] = lcv_u;
                Elg[q] = gv;
            }
        }
        
        GElements[glob_id]=El;
        
        if(El.size()==4)
        {
            Tetras[glob_id]=El;
            GTetras[glob_id]=Elg;
        }
        if(El.size()==6)
        {
            Prisms[glob_id]=El;
            GPrisms[glob_id]=Elg;
        }
        if(El.size()==8)
        {
            Hexes[glob_id]=El;
            GHexes[glob_id]=Elg;
        }
        El.clear();
        Elg.clear();
        elid++;
    }
    
    pDom->Elements        = Elements;
    pDom->Hexes           = Hexes;
    pDom->Prisms          = Prisms;
    pDom->Tetras          = Tetras;
    pDom->GHexes           = GHexes;
    pDom->GPrisms          = GPrisms;
    pDom->GTetras          = GTetras;
    pDom->LocElem2LocNode = locelem2locnode;
    pDom->loc_part_verts  = loc_part_verts;
    pDom->glob_part_verts = glob_part_verts;
    pDom->gv2lpv          = gv2lpv;
    pDom->lv2gpv          = lv2gpv;
    pDom->vert2elem       = vert2elem;
    pDom->gv2lpartv       = gv2lpartv;
    pDom->lpartv2gv       = lpartv2gv;
    
}


std::map<int,std::map<int,double> > Partition::getNode2NodeMap()
{
    std::map<int,std::map<int,double> > node2node;
    std::map<int,std::vector<int> >::iterator ieet;
    int gvidt,gvid,gel,Nv,Nf;
    Vert* V0 = new Vert;
    Vert* V1 = new Vert;
    
    
    for(ieet=iee_part_map->i_map.begin();ieet!=iee_part_map->i_map.end();ieet++)
    {
        Nv  = LocElem2Nv[ieet->first];
        Nf  = LocElem2Nf[ieet->first];
        
        for(int i=0;i<Nv;i++)
        {
            gvid = ien_part_map->i_map[ieet->first][i];
            
            for(int j=0;j<Nv;j++)
            {
                gvidt = ien_part_map->i_map[ieet->first][j];
                if(gvid!=gvidt && node2node[gvid].find(gvidt)==node2node[gvid].end())
                {
                    
                    int lvid  = GlobalVert2LocalVert[gvid];
                    int lvidt = GlobalVert2LocalVert[gvidt];

                    V0->x = LocalVerts[lvid]->x;
                    V0->y = LocalVerts[lvid]->y;
                    V0->z = LocalVerts[lvid]->z;
                    
                    V1->x = LocalVerts[lvidt]->x;
                    V1->y = LocalVerts[lvidt]->y;
                    V1->z = LocalVerts[lvidt]->z;
                    
//                    double dist = ComputeEdgeLength(V0,V1);
                    double dist = sqrt((V0->x-V1->x)*(V0->x-V1->x)
                                +(V0->y-V1->y)*(V0->y-V1->y)
                                +(V0->z-V1->z)*(V0->z-V1->z));
                    node2node[gvid].insert(std::pair<int,double>(gvidt,dist));
                    
                    
                }
            }
            
            
//            for(int k=0;k<Nf;k++)
//            {
//                gel = iee_part_map->i_map[ieet->first][k];
//
//                if(gel<NelGlob)
//                {
//                    int Nvv  = LocElem2Nv[gel];
//
//                    for(int j=0;j<Nvv;j++)
//                    {
//                        gvidt = ien_part_map->i_map[gel][j];
//                        if(gvidt!=gvid && node2node[gvid].find(gvidt)==node2node[gvid].end())
//                        {
//                            int lvid  = GlobalVert2LocalVert[gvid];
//                            int lvidt = GlobalVert2LocalVert[gvidt];
//
//                            V0->x = LocalVerts[lvid]->x;
//                            V0->y = LocalVerts[lvid]->y;
//                            V0->z = LocalVerts[lvid]->z;
//
//                            V1->x = LocalVerts[lvidt]->x;
//                            V1->y = LocalVerts[lvidt]->y;
//                            V1->z = LocalVerts[lvidt]->z;
//
//                            //double dist = ComputeEdgeLength(V0,V1);
//                            double dist = sqrt((V0->x-V1->x)*(V0->x-V1->x)
//                                        +(V0->y-V1->y)*(V0->y-V1->y)
//                                        +(V0->z-V1->z)*(V0->z-V1->z));
//
//                            node2node[gvid].insert(std::pair<int,double>(gvidt,dist));
//                        }
//                    }
//                }
//            }
        }
    }
    delete V0;
    delete V1;
    
    return node2node;
}

std::map<int,double> Partition::ReduceFieldToAllVertices(std::map<int,double> UaddAdj)
{
    std::map<int,double> Uvm;
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    //std::map<int,std::vector<int> > gV2gE = globElem2globVerts;
    
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        double sum = 0.0;
        double avg = 0.0;
        
        int gv = itm->first;
        
        for(int q=0;q<globVerts2globElem[gv].size();q++)
        {
            int gEl = globVerts2globElem[gv][q];
            sum     = sum + UaddAdj[gEl];
        }
        avg = sum/globVerts2globElem[gv].size();
        Uvm[gv]=avg;

        im++;
    }
    
    return Uvm;

}





std::map<int,Array<double>* > Partition::ReduceStateVecToAllVertices(std::map<int,Array<double>* > UaddAdj, int nvar)
{
    std::map<int,Array<double>* > Uvm;
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    
    double* sum = new double[nvar];
//    Array<double>* avg = new Array<double>(nvar,1);
    
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        int gv = itm->first;
        
        for(int l=0;l<nvar;l++)
        {
            sum[l] = 0.0;
        }
        
        for(int q=0;q<globVerts2globElem[gv].size();q++)
        {
            int gEl = globVerts2globElem[gv][q];
            for(int l=0;l<nvar;l++)
            {
                sum[l] = sum[l] + UaddAdj[gEl]->getVal(l,0);
            }
        }
        Array<double>* avg = new Array<double>(nvar,1);
        for(int l=0;l<nvar;l++)
        {
            avg->setVal(l,0,sum[l]/globVerts2globElem[gv].size());
        }
        
        Uvm[gv]=avg;

        im++;
    }
    
    
    delete[] sum;
    
    return Uvm;

}






std::map<int,double> Partition::ReduceFieldToVertices(std::map<int,double> Uelem)
{
    std::vector<double> Uv;
    std::map<int,double> Uvm;
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        double sum = 0.0;
        double avg = 0.0;
        
        int gv = itm->first;
        
        for(int q=0;q<itm->second.size();q++)
        {
            int gEl = itm->second[q];
            sum = sum + Uelem[gEl];
        }
        
        avg = sum/itm->second.size();
        Uv.push_back(avg);
        Uvm[gv]=avg;
        //std::cout <<"hhh "<< gv << " " << gv2  << " " << avg<<std::endl;

        im++;
    }

    return Uvm;

}


std::map<int,Array<double>* > Partition::ReduceMetricToVertices(std::map<int,Array<double>* > Telem)
{
    std::map<int,Array<double>*> Tmet_v;
    
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    double sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
    double sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
    double sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
        sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
        sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
        double avg = 0.0;
        int gv = itm->first;
        Array<double>* Tv = new Array<double>(3,3);
        for(int q=0;q<itm->second.size();q++)
        {
            int gEl = itm->second[q];
            
            sum00 = sum00 + Telem[gEl]->getVal(0,0);
            sum01 = sum01 + Telem[gEl]->getVal(0,1);
            sum02 = sum02 + Telem[gEl]->getVal(0,2);
            
            sum10 = sum10 + Telem[gEl]->getVal(1,0);
            sum11 = sum11 + Telem[gEl]->getVal(1,1);
            sum12 = sum12 + Telem[gEl]->getVal(1,2);
            
            sum20 = sum20 + Telem[gEl]->getVal(2,0);
            sum21 = sum21 + Telem[gEl]->getVal(2,1);
            sum22 = sum22 + Telem[gEl]->getVal(2,2);
        }
        
        Tv->setVal(0,0,sum00/itm->second.size());
        Tv->setVal(0,1,sum01/itm->second.size());
        Tv->setVal(0,2,sum02/itm->second.size());
        
        Tv->setVal(1,0,sum10/itm->second.size());
        Tv->setVal(1,1,sum11/itm->second.size());
        Tv->setVal(1,2,sum12/itm->second.size());
        
        Tv->setVal(2,0,sum20/itm->second.size());
        Tv->setVal(2,1,sum21/itm->second.size());
        Tv->setVal(2,2,sum22/itm->second.size());
        
        Tmet_v[gv]=Tv;

        im++;
    }
    
    return Tmet_v;
}


//std::map<int,Array<double>*> Partition::ReduceMetricToAllVertices(std::map<int,Array<double>* > Telem)
//{
//    std::vector<double> Uv;
//    std::map<int,Array<double>*> Tmet_v;
//    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
//    std::vector<int> glob_part_verts = pDom->glob_part_verts;
//    std::map<int,int> lv2gpv = pDom->lv2gpv;
//    int im = 0;
//    int tel=0;
//    std::map<int,std::vector<int> >::iterator itm;
//    double sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
//    double sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
//    double sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
//    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
//    {
//        double sum = 0.0;
//        double avg = 0.0;
//        
//        sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
//        sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
//        sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
//        
//        int gv = itm->first;
//        Array<double>* Tv = new Array<double>(3,3);
//
//        for(int q=0;q<globVerts2globElem[gv].size();q++)
//        {
//            int gEl = globVerts2globElem[gv][q];
//            sum00 = sum00 + Telem[gEl]->getVal(0,0);
//            sum01 = sum01 + Telem[gEl]->getVal(0,1);
//            sum02 = sum02 + Telem[gEl]->getVal(0,2);
//            
//            sum10 = sum10 + Telem[gEl]->getVal(1,0);
//            sum11 = sum11 + Telem[gEl]->getVal(1,1);
//            sum12 = sum12 + Telem[gEl]->getVal(1,2);
//            
//            sum20 = sum20 + Telem[gEl]->getVal(2,0);
//            sum21 = sum21 + Telem[gEl]->getVal(2,1);
//            sum22 = sum22 + Telem[gEl]->getVal(2,2);
//        }
//        avg = sum/globVerts2globElem[gv].size();
//        
//        Tv->setVal(0,0,sum00/itm->second.size());
//        Tv->setVal(0,1,sum01/itm->second.size());
//        Tv->setVal(0,2,sum02/itm->second.size());
//        
//        Tv->setVal(1,0,sum10/itm->second.size());
//        Tv->setVal(1,1,sum11/itm->second.size());
//        Tv->setVal(1,2,sum12/itm->second.size());
//        
//        Tv->setVal(2,0,sum20/itm->second.size());
//        Tv->setVal(2,1,sum21/itm->second.size());
//        Tv->setVal(2,2,sum22/itm->second.size());
//        
//        //Uv.push_back(avg);
//        Tmet_v[gv]=Tv;
//        im++;
//    }
//    
//    return Tmet_v;
//
//}


Domain* Partition::getPartitionDomain()
{
    return pDom;
}

std::vector<int> Partition::getLocElem()
{
    return Loc_Elem;
}
std::vector<int> Partition::getLocElemNv()
{
    return Loc_Elem_Nv;
}
std::map<int,int> Partition::getLocElem2Nv()
{
    return LocElem2Nv;
}
std::map<int,int> Partition::getLocElem2Nf()
{
    return LocElem2Nf;
}
std::vector<double> Partition::getLocElemVaria()
{
    return Loc_Elem_Varia;
}
std::vector<int> Partition::getLocAndAdjElem()
{
    return LocAndAdj_Elem;
}

int Partition::getnLoc_Elem()
{
    return nLoc_Elem;
}

int Partition::getnLoc_Verts()
{
    return nLoc_Verts;
}

int* Partition::getXadj()
{
    return xadj;
}
int* Partition::getAdjcny()
{
    return adjcny;
}
ParArray<int>* Partition::getLocalPartition()
{
    return part;
}
Array<int>* Partition::getGlobalPartition()
{
    return part_global;
}
std::vector<Vert*> Partition::getLocalVerts()
{
    return LocalVerts;
}
Vert* Partition::getLocalVert(int v_loc_id)
{
    return LocalVerts[v_loc_id];
}
std::vector<std::vector<int> > Partition::getLocalElem2GlobalVert()
{
    return LocalElem2GlobalVert;
}
std::vector<std::vector<int> > Partition::getLocalElem2LocalVert()
{
    return LocalElem2LocalVert;
}
std::map<int,int> Partition::getLocalVert2GlobalVert()
{
    return LocalVert2GlobalVert;
}
std::map<int,int> Partition::getGlobalVert2LocalVert()
{
    return GlobalVert2LocalVert;
}
std::map<int,int> Partition::getLocalFace2GlobalFace()
{
    return LocalFace2GlobalFace;
}
std::map<int,int> Partition::getGlobalFace2LocalFace()
{
    return GlobalFace2LocalFace;
}
//std::map<int, std::vector<int> > Partition::getglobElem2localFaces()
//{
//    return globElem2localFaces;
//}
std::map<int, std::vector<int> > Partition::getglobElem2globFaces()
{
    return globElem2globFaces;
}
std::map<int, std::vector<int> > Partition::getglobFace2GlobalElements()
{
    return globFace2GlobalElements;
}
std::set<int> Partition::getElemSet()
{
    return elem_set;
}
std::set<int> Partition::getLocElemSet()
{
    return loc_r_elem_set;
}
//std::vector<double> Partition::getUelem()
//{
//    return U0Elem;
//}
//double Partition::getU0atGlobalElem(int gelem)
//{
//    int elem = GlobalElement2LocalElement[gelem];
//    return U0Elem[elem];
//}
Array<double>* Partition::getUvert()
{
    return U0Vert;
}
std::map<int,std::vector<int> > Partition::getGlobElem2LocVerts()
{
    return globElem2locVerts;
}
std::map<int,std::vector<int> > Partition::getGlobElem2GlobVerts()
{
    return globElem2globVerts;
}
std::map<int,std::vector<int> > Partition::getGlobVert2GlobElem()
{
    return globElem2globVerts;
}
std::map<int,int> Partition::getGlobalElement2LocalElement()
{
    return GlobalElement2LocalElement;
}
std::map<int,int> Partition::getLocalElement2GlobalElement()
{
    return LocalElement2GlobalElement;
}
ParallelState* Partition::getParallelState()
{
    return ien_pstate;
}
i_part_map* Partition::getIEEpartmap()
{
    return iee_part_map;
}
i_part_map* Partition::getIEFpartmap()
{
    return ief_part_map;
}
i_part_map* Partition::getIENpartmap()
{
    return ien_part_map;
}
i_part_map* Partition::getIFNpartmap()
{
    return ifn_part_map;
}
i_part_map* Partition::getIF_Nvpartmap()
{
    return if_Nv_part_map;
}
i_part_map* Partition::getIFEpartmap()
{
    return ife_part_map;
}
i_part_map* Partition::getIFREFpartmap()
{
    return if_ref_part_map;
}
ParallelState* Partition::getXcnParallelState()
{
    return xcn_pstate;
}
ParallelState* Partition::getIenParallelState()
{
    return ien_pstate;
}
