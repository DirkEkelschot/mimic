#include "adapt_partition.h"

Partition::Partition(ParArray<int>* ien, ParArray<int>* iee, ParArray<int>* ief, ParArray<int>* ifn, ParArray<int>* ife, ParArray<int>* if_ref,  ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, ParallelState* ife_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    ien_pstate = ien_parstate;
    xcn_pstate = xcn_parstate;
    // This function computes the xadj and adjcny array and the part array which determines which element at current rank should be sent to other ranks.
    
    double t0 = MPI_Wtime();
    DeterminePartitionLayout(ien, pstate_parmetis, ien_pstate, comm);
    double t1 = MPI_Wtime();
    double time_layout = t1-t0;
    double max_time_layout = 0.0;
    MPI_Allreduce(&time_layout, &max_time_layout, 1, MPI_DOUBLE, MPI_MAX, comm);
    
    eloc = 0;
    vloc = 0;
    floc = 0;

    // This function takes care of the send and receive operations in order to send the appropriate elements and corresponding vertices to the appropriate rank.
    // These operations are based on the fact that part holds the desired partitioning of the elements. the spread of the vertices is based on the fact that all the vertices stored in xcn are distributed "uniformly";

    DetermineElement2ProcMap(ien, ief, ien_pstate, xcn, xcn_parstate, U, comm);

    //DetermineAdjacentElement2ProcMap(ien, ief, part, ien_pstate, xcn, xcn_parstate, U, comm);
    iee_part_map    = getElement2EntityPerPartition(iee,    ien_pstate,   comm);
    ief_part_map    = getElement2EntityPerPartition(ief,    ien_pstate,   comm);
    ien_part_map    = getElement2EntityPerPartition(ien,    ien_pstate,   comm);
    
    ifn_part_map    = getFace2EntityPerPartition(ifn    ,    ife_parstate,   comm);
    ife_part_map    = getFace2EntityPerPartition(ife    ,    ife_parstate,   comm);
    if_ref_part_map = getFace2EntityPerPartition(if_ref ,    ife_parstate,   comm);
    

    
    DetermineAdjacentElement2ProcMapUS3D(ien, iee_part_map->i_map, part, ien_pstate, xcn, xcn_parstate, U, comm);

    nLocAndAdj_Elem = LocAndAdj_Elem.size();
    
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void Partition::DeterminePartitionLayout(ParArray<int>* ien, ParallelState_Parmetis* pstate_parmetis, ParallelState* ien_parstate, MPI_Comm comm)
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
    
    /*
    for(int u=0;u<nloc;u++)
    {
        part_arr[u] = rank;
    }
    */

    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(),
                           xadj_par, adjncy_par,
                               elmwgt, adjwgt,
                       vsize, wgtflag,
                   numflag, ncon, nparts,
                   tpwgts, ubvec, itr, options,
                   &edgecut, part_arr, &comm);
    
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
                   ien_pstate->getNlocs(),
                   ien_pstate->getOffsets(),
                   MPI_INT,comm);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void Partition::DetermineElement2ProcMap(ParArray<int>* ien, ParArray<int>* ief, ParallelState* ien_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm)
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
    
    int ien_o = part->getOffset(rank);
    double varia = 0.0;
    int not_on_rank=0;
    int on_rank = 0;
    std::set<int> requested_elements;
    int* new_offsets = new int[size];
    
    for(i=0;i<size;i++)
    {
        new_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }
    
    int* new_offsets2 = new int[size];
    for(int i=0;i<size;i++)
    {
       new_offsets2[i] = ien_pstate->getOffsets()[i]-1;
    }
    for(i=0;i<part->getNrow();i++)
    {
        p_id  = part->getVal(i,0);
        el_id = part->getOffset(rank)+i;
        varia   = U->getVal(i,0);
        //std::cout << "rho = " << rho << std::endl;
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            varia_to_send_to_ranks[p_id].push_back(varia);
            
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
                        vloc_tmp++;
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
                        vloc_tmp++;
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
                        floc_tmp++;
                    }
                }
            }
            loc_elem.push_back(el_id);
            loc_varia.push_back(varia);
            elem_set.insert(el_id);
            loc_elem_set.insert(el_id);
            
            on_rank++;
        }
    }
    
    
    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);
    
    std::map<int,std::vector<int> >  part_tot_recv_elIDs_map;
    std::map<int,std::vector<double> >  part_tot_recv_varias_map;
    
    std::map<int,std::vector<int> >  TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_f_map;
    std::map<int,std::vector<int> >::iterator it;
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
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&varia_to_send_to_ranks[it->first][0], n_req, MPI_DOUBLE, dest, 20000+100+dest*2, comm);
                
                //MPI_Send(&n_req_v, 1, MPI_INT, dest, 9000+dest, comm);
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<double> part_recv_varia(n_req_recv);
            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_vrt_id(n_req_recv*8);
            std::vector<int>    part_recv_face_id(n_req_recv*6);
            
            
            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_varia[0], n_req_recv, MPI_DOUBLE, q, 20000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            MPI_Recv(&part_recv_vrt_id[0], n_req_recv*8, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv*6, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            
            TotRecvElement_IDs_v_map[q] = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q] = part_recv_face_id;
            part_tot_recv_elIDs_map[q]  = part_recv_el_id;
            part_tot_recv_varias_map[q] = part_recv_varia;
            
        }
    }
    
    std::vector<int> TotRecvElement_IDs;
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
        
        for(int k=0;k<8;k++)
        {
            int v_id_n = TotRecvVerts_IDs[cnt_v];
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
                    vloc_tmp++;
                }
            }
            cnt_v++;
            
            if(k<6)
            {
                int f_id_n = TotRecvFaces_IDs[cnt_f];
                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
                {
                    unique_faceIDs_on_rank_set.insert(f_id);
                    faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                    floc_tmp++;
                }
                cnt_f++;

            }
            
        }
        //part_elem2verts.push_back(elem);
        elem_set.insert(TotRecvElement_IDs[i]);
        loc_elem_set.insert(TotRecvElement_IDs[i]);
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
                offset_xcn        = xcn_parstate->getOffset(rank);
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
    double* part_verts_arr = new double[3*(vloc_tmp+vfor)];

    int gvid=0;
    int lvid=0;
    
    for(m=0;m<vloc_tmp;m++)
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
    
    int o = 3*vloc_tmp;
    m = 0;
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
    //NlocElem             = loc_elem.size()+Nel_extra+itel;
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
    
    for(m=0;m<loc_elem.size();m++)
    {
        el_id   = loc_elem[m];
        varia_v = loc_varia[m];
        
        Loc_Elem.push_back(el_id);
        LocAndAdj_Elem.push_back(el_id);
        Loc_Elem_Varia.push_back(varia_v);
        LocAndAdj_Elem_Varia.push_back(varia_v);
        //U0Elem->setVal(m,0,rho_v);
        //ElemPart->setVal(m,0,el_id);
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        
        for(int p=0;p<8;p++)
        {
            //GlobalFace2LocalFace
            //LocalVert2GlobalVert
            glob_v = ien->getVal(el_id-ien_o,p);
            loc_v  = GlobalVert2LocalVert[glob_v];
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            //LocalElem2GlobalVert->setVal(m,p,glob_v);
            //LocalElem2LocalVert->setVal(m,p,loc_v);
            //collect_var[loc_v].push_back(rho_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            
            if(p<6)
            {
                glob_f = ief->getVal(el_id-ien_o,p);
                loc_f  = GlobalFace2LocalFace[glob_f];
                globElem2localFaces[el_id].push_back(loc_f);
                globElem2globFaces[el_id].push_back(glob_f);
                globFace2GlobalElements[glob_f].push_back(el_id);
            }
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
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        varia_v = TotRecvElement_varia[m];
        Loc_Elem.push_back(el_id);
        LocAndAdj_Elem.push_back(el_id);
        Loc_Elem_Varia.push_back(varia_v);
        LocAndAdj_Elem_Varia.push_back(varia_v);

        for(int p=0;p<8;p++)
        {
            glob_v = TotRecvVerts_IDs[cnv];
            loc_v = GlobalVert2LocalVert[glob_v];
            //LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            //LocalElem2LocalVert->setVal(m+o,p,loc_v);
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            //collect_var[loc_v].push_back(rho_v);
            cnv++;
            if(p<6)
            {
                glob_f = TotRecvFaces_IDs[cnf];
                loc_f  = GlobalFace2LocalFace[glob_f];
                globElem2localFaces[el_id].push_back(loc_f);
                globElem2globFaces[el_id].push_back(glob_f);
                globFace2GlobalElements[glob_f].push_back(el_id);
                cnf++;
            }
        }
        
        LocalElem2GlobalVert.push_back(tmp_globv);
        LocalElem2LocalVert.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }
    
    nLoc_Elem = Loc_Elem.size();
    vloc = LocalVerts.size();
    floc = cnf;
    delete[] new_offsets;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Partition::DetermineAdjacentElement2ProcMapUS3D(ParArray<int>* ien, std::map<int,std::vector<int> > iee_vec, ParArray<int>* part, ParallelState* ien_parstate, ParArray<double>* xcn, ParallelState* xcn_parstate, Array<double>* U, MPI_Comm comm)
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
    
    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_vert;
    int* new_offsets = new int[size];
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = xcn_parstate->getOffsets()[i]-1;
    }
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    int Nel = part_global->getNrow();
    for(int i=0;i<Loc_Elem.size();i++)
    {
        int elId = Loc_Elem[i];

        int k     = 0;
        for(int j=0;j<6;j++)
        {
            int adjEl_id = iee_vec[elId][j];
            
            if((elem_set.find(adjEl_id)==elem_set.end()) && adjEl_id<Nel)
            {
                p_id = part_global->getVal(adjEl_id,0);
                if(p_id != rank)
                {
                    adj_elements[p_id].push_back(adjEl_id);
                    elem_set.insert(adjEl_id);
                    req_elem[p_id].push_back(adjEl_id);
                    itel++;
                }
            }
        }
    }
    
    adj_schedule = DoScheduling(req_elem,comm);
    
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
            TotAdj_El_IDs.push_back(adj_id);
            //TotAdj_Rhos.push_back(U->getVal(adj_id-offset_new,0));
            //send_adj_rhos[dest].push_back(U->getVal(adj_id-offset_new,0));
            //send_adj_rhos[dest].push_back(U0Elem[lelem]);
            for(int k=0;k<8;k++)
            {
                v_id = globElem2globVerts[adj_id][k];
                //v_id = ien->getVal(adj_id-offset_new,k);
                send_adj_verts_IDs[dest].push_back(v_id);
            }
        
            for(int k=0;k<6;k++)
            {
                f_id = globElem2globFaces[adj_id][k];
                send_adj_faces_IDs[dest].push_back(f_id);
            }
        }

        TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn = 0;
    std::map<int,int  > recv_adj_back_Nverts;
    std::map<int,int* > recv_adj_back_verts_ids;
    std::map<int,int  > recv_adj_back_Nfaces;
    std::map<int,int* > recv_adj_back_faces_ids;
    //std::map<int,int > recv_adj_back_Nrhos;
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

            //recv_adj_back_Nrhos[q]      = n_adj_rho_recv_back;
            //recv_adj_back_rhos[q]       = recv_adj_back_arr_rho;

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
    std::map<int,std::vector<int> >::iterator itm_el;
    for(itm_el=adj_elements.begin();itm_el!=adj_elements.end();itm_el++)
    {
        //TotNrho_adj_recv = TotNrho_adj_recv+itm->second;
        for(int i=0;i<itm_el->second.size();i++)
        {
            //adj_rhos.push_back(recv_adj_back_rhos[itm->first][i]);
            
            adj_elements_vec.push_back(adj_elements[itm_el->first][i]);
            
        }
    }

    //std::cout << "TotNelem_adj_recv " << TotNelem_adj_recv << " should be equal to " << TotNrho_adj_recv << std::endl;
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
                    vloc_tmp++;
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
                    floc_tmp++;
                }
                cnt_f_adj++;
            }
        }
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
                   offset_xcn        = xcn_parstate->getOffset(rank);
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
       double* part_verts_arr = new double[3*(vloc_tmp+vfor)];


       int gvid=0;
       int lvid=vloc;
       int m=0;
       for(m=0;m<vloc_tmp;m++)
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
   
       int o = 3*vloc_tmp;
       m = 0;
       int u = 0;
       for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
       {
           int Nv = recv_back_Nverts[it_f->first];
          
           for(u=0;u<Nv;u++)
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
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        //rho_v = adj_rhos[m];

        //U0Elem.push_back(rho_v);
        LocAndAdj_Elem.push_back(el_id);
        //LocAndAdj_Elem_Varia.push_back(varia_v);
//      U0Elem->setVal(m+o,0,rho_v);
//      ElemPart->setVal(m+o,0,el_id);

        for(int p=0;p<8;p++)
        {
            glob_v = adj_verts[cnv];
            loc_v  = GlobalVert2LocalVert[glob_v];
            //LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            //LocalElem2LocalVert->setVal(m+o,p,loc_v);
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            //collect_var[loc_v].push_back(rho_v);
            cnv++;
            if(p<6)
            {
                glob_f = adj_faces[cnf];
                loc_f  = GlobalFace2LocalFace[glob_f];
                globElem2localFaces[el_id].push_back(loc_f);
                globElem2globFaces[el_id].push_back(glob_f);
                globFace2GlobalElements[glob_f].push_back(el_id);
                cnf++;
            }
        }
        LocalElem2GlobalVert.push_back(tmp_globv);
        LocalElem2LocalVert.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }
     
    //nLoc_Elem = U0Elem.size();
    
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
    
    
    delete[] part_verts_arr;
    delete[] new_offsets;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


std::map<int,double> Partition::CommunicateAdjacentDataUS3D(std::map<int,double> U, MPI_Comm comm)
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
    
    int ien_o = ien_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Elems;
    int* new_offsets = new int[size];
    std::map<int,double> U_loc;
    
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
        else
        {
            U_loc[el_req]=U[el_req];
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
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_ghost_IDs;
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
                int offset_iee          = ien_pstate->getOffset(rank);
                
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
             
            double* recv_back_iee_arr   = new double[n_recv_back*6];
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
            U_loc[el_id] = recv_back_iee[iter->first][s];
        }
    }
    
    delete[] new_offsets;
    
    return U_loc;
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



i_part_map* Partition::getElement2EntityPerPartition(ParArray<int>* iee, ParallelState* ien_pstate, MPI_Comm comm)
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
    int* new_offsets = new int[size];
    std::map<int,std::vector<int> > iee_loc;
    std::map<int,std::vector<int> > iee_loc_inv;
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ien_pstate->getOffsets()[i]-1;
    }
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    int ncol = iee->getNcol();
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    for(int i=0;i<Loc_Elem.size();i++)
    {
        int el_req = Loc_Elem[i];
        
        r = FindRank(new_offsets,size,el_req);
        
        if(r != rank)
        {
            rank2req_Elems[r].push_back(el_req);
        }
        else
        {
            for(int j=0;j<ncol;j++)
            {
                iee_loc[el_req].push_back(iee->getVal(el_req-ien_o,j));
                iee_loc_inv[iee->getVal(el_req-ien_o,j)].push_back(el_req);
            }
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

                
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_ghost_IDs;
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
                double* iee_send        = new double[ne_send*ncol];
                int offset_iee          = ien_pstate->getOffset(rank);
                
                for(int u=0;u<it->second.size();u++)
                {
                    for(int h=0;h<ncol;h++)
                    {
                        iee_send[u*ncol+h] = iee->getVal(it->second[u]-offset_iee,h);
                    }
                }

                int dest = it->first;
                MPI_Send(&ne_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], ne_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&iee_send[0], ne_send*ncol, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                delete[] iee_send;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_iee_arr   = new double[n_recv_back*ncol];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_iee_arr[0], n_recv_back*ncol, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

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
            for(int r=0;r<ncol;r++)
            {
                iee_loc[el_id].push_back(recv_back_iee[iter->first][s*ncol+r]);
                iee_loc_inv[recv_back_iee[iter->first][s*ncol+r]].push_back(el_id);
                
            }
        }
        ntotal=ntotal+L;
    }
    
    delete[] new_offsets;
    
    iee_p_map->i_map = iee_loc;
    iee_p_map->i_inv_map = iee_loc_inv;
    
    return iee_p_map;
}




i_part_map* Partition::getFace2EntityPerPartition(ParArray<int>* ife, ParallelState* ife_pstate, MPI_Comm comm)
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
    //std::cout << xcn->getOffset(rank) << " " << xcn_pstate->getOffset(rank) << std::endl;
    
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
    std::map<int,std::vector<int> > send_ghost_IDs;
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
    
//    for(int i=0;i<ife_loc.size();i++)
//    {
//        for(int j=0;j<ife_loc[0].size();j++)
//        {
//            std::cout << ife_loc[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
    delete[] new_offsets;
    
    ife_p_map->i_map = ife_loc;
    ife_p_map->i_inv_map = ife_loc_inv;
    
    return ife_p_map;
}






std::vector<int> Partition::getLocElem()
{
    return Loc_Elem;
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
std::vector<Vert> Partition::getLocalVerts()
{
    return LocalVerts;
}
Vert Partition::getLocalVert(int v_loc_id)
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
    return loc_elem_set;
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
