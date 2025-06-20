


#include "adapt_partobject_lite.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"
#include "adapt_elements.h"

PartObjectLite::PartObjectLite(std::map<int,std::vector<int> > Elem2Vert_uniform,
                               std::map<int,std::vector<double> > vertices_i,
                               std::map<int,int> eltype_map,
                               std::vector<int> elTypes,
                               FaceSetPointer allbcFaces,
                               int nE,
                               int nV,
                               MPI_Comm comm)
                               
{
    MPI_Comm_size(comm, &m_size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &m_rank);

    m_Ne_glob = nE;
    m_Nv_glob = nV;
    //std::cout << "Elem2Vert_uniform " << Elem2Vert_uniform.size() << std::endl;
    DeterminePartitionLayout(Elem2Vert_uniform,elTypes,comm);
    
    DetermineElement2ProcMap(Elem2Vert_uniform, 
                             vertices_i,
                             eltype_map, 
                             comm);

    //std::map<int,std::vector<int> > Rank2Elem = CommunicateAdjacencyInfoLocalPartition(comm);

    ComputeFaceMap(comm, allbcFaces);

    BuildPMMGCommunicationData(comm, allbcFaces);

}







void PartObjectLite::DeterminePartitionLayout(std::map<int,std::vector<int> > Elem2Vert_uniform,
                                            std::vector<int> elTypes,
                                            MPI_Comm comm)
{
    int root = 0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int nrow    = Elem2Vert_uniform.size();
    int nvpEL   = Elem2Vert_uniform.begin()->second.size();
    int nloc    = nrow;

    //std::cout << " elTypes " << elTypes[0] << " " << elTypes[1] << " " << elTypes[2] << " " << nvpEL << std::endl;

    ParallelState_Parmetis_Lite* pstate_parmetis = new ParallelState_Parmetis_Lite(Elem2Vert_uniform,  elTypes, comm);
    
    //=================================================================
    //=================================================================
    //=================================================================

    //ParallelState_Parmetis* pstate_parmetis2 = new ParallelState_Parmetis(ien,comm,8);
    //
    idx_t numflag_[]        = {0};
    idx_t *numflag          = numflag_;
    idx_t ncommonnodes_[]   = {pstate_parmetis->getNcommonNodes()};
    idx_t *ncommonnodes     = ncommonnodes_;
    int edgecut             = 0;
    idx_t *xadj_par         = NULL;
    idx_t *adjncy_par       = NULL;
    idx_t options_[]        = {0, 0, 0};
    idx_t *options          = options_;
    idx_t wgtflag_[]        = {2};
    idx_t *wgtflag          = wgtflag_;
    real_t ubvec_[]         = {1.1};
    real_t *ubvec           = ubvec_;

    std::vector<int> elmwgt = pstate_parmetis->getElmWgt();

    int np                  = size;
    idx_t ncon_[]           = {1};
    idx_t *ncon             = ncon_;
    real_t *tpwgts          = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
    tpwgts[i] = 1.0/np;
    }

    idx_t nparts_[]         = {np};
    idx_t *nparts           = nparts_;
    std::vector<int> part_arr(nloc,0);
    real_t itr_[]           = {1.05};
    real_t *itr             = itr_;

    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;

    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist().data(),
    pstate_parmetis->getEptr().data(),
    pstate_parmetis->getEind().data(),
    numflag,ncommonnodes,
    &xadj_par,&adjncy_par,&comm);

    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist().data(),
    xadj_par,
    adjncy_par,
    pstate_parmetis->getElmWgt().data(), NULL, wgtflag, numflag,
    ncon, nparts,
    tpwgts, ubvec, options,
    &edgecut, part_arr.data(), &comm);

    m_part = std::vector<int>(Elem2Vert_uniform.size(),0);
    int nElemTotal = pstate_parmetis->getNtotalElem();


    std::vector<int> part_arr_ell_id(nloc,0);
    std::vector<int> part_arr_ell_Rid(nloc,0);

    std::map<int,std::vector<int> >::iterator itmiv;
    int i = 0;

    for(itmiv=Elem2Vert_uniform.begin();itmiv!=Elem2Vert_uniform.end();itmiv++)
    {
    int gid                 = itmiv->first;
    part_arr_ell_id[i]      = gid;
    part_arr_ell_Rid[i]     = part_arr[i];
    m_partMap[gid]          = part_arr[i];
    i++;
    }

    std::vector<int> m_partGlobalRoot_vec;
    std::vector<int> m_partGlobalRoot_Rid_vec;
    if ( rank == root) 
    {
    m_partGlobalRoot_vec.resize(nElemTotal,0);
    m_partGlobalRoot_Rid_vec.resize(nElemTotal,0);
    } 

    MPI_Gatherv(&part_arr_ell_id.data()[0],
    nloc, MPI_INT,
    &m_partGlobalRoot_vec.data()[0],
    pstate_parmetis->getNlocs().data(),
    pstate_parmetis->getElmdist().data(),
    MPI_INT,root,comm);

    MPI_Gatherv(&part_arr_ell_Rid.data()[0],
    nloc, MPI_INT,
    &m_partGlobalRoot_Rid_vec.data()[0],
    pstate_parmetis->getNlocs().data(),
    pstate_parmetis->getElmdist().data(),
    MPI_INT,root,comm);

    for(int i=0;i<m_partGlobalRoot_vec.size();i++)
    {
    int eid = m_partGlobalRoot_vec[i];
    int rid = m_partGlobalRoot_Rid_vec[i];

    m_partGlobalRoot[eid] = rid;
    }

    m_nElemGlobalPart = m_partGlobalRoot.size();

    MPI_Bcast(&m_nElemGlobalPart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] xadj_par;
    delete[] adjncy_par;
    /**/
}



void PartObjectLite::DetermineElement2ProcMap(std::map<int,std::vector<int> >   Elem2Vert_uniform, 
    std::map<int,std::vector<double> > vertices_i,
    std::map<int,int>                  Elem2Type_uniform,
    MPI_Comm comm)
{   

    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > nvPerElms_to_send_to_ranks;
    std::map<int,std::vector<int> > ndataPerElms_to_send_to_ranks;
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
    std::vector<std::vector<int> >  part_elem2verts;
    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > elmtype_to_send_to_ranks;
    ParallelState* ife_pstate  = new ParallelState(m_Nf_glob,comm);
    ParallelState* vert_pstate = new ParallelState(m_Nv_glob,comm);

    std::map<int,std::vector<int> > rank2req_vert;

    std::vector<int> vertIDs_on_rank;
    std::vector<int> part_v;
    
    int r     = 0;
    int lv_id = 0;
    int lf_id = 0;
    int f_id  = 0;
    int ea_id = 0;

    int not_on_rank=0;
    int on_rank = 0;
    int* new_V_offsets = new int[size];
    for(i=0;i<size;i++)
    {
        new_V_offsets[i] = vert_pstate->getOffsets()[i]-1;
    }
    
    int nvPerEl;
    int nfPerEl;
    int ndataPerEl;
    int tett = 0;
    std::map<int,int>::iterator itmii;
    
    std::vector<double> loc_r_data;
    std::vector<int> loc_r_Ndata;
    std::vector<int> loc_r_elem;
    std::vector<int> loc_r_nf_elem;
    std::vector<int> loc_r_nv_elem;

    for(itmii=m_partMap.begin();itmii!=m_partMap.end();itmii++)
    {
        el_id       = itmii->first;
        p_id        = itmii->second;
        int el_type = Elem2Type_uniform[el_id];
        nvPerEl     = Elem2Vert_uniform[el_id].size();
        // if(el_type != 10)
        // {
        //     std::cout << "el_id " << el_id << " p_id " << p_id << " nvPerEl " << nvPerEl << " el_type " << el_type << " " << nvPerEl << std::endl;
        // }
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            nvPerElms_to_send_to_ranks[p_id].push_back(nvPerEl);
            elmtype_to_send_to_ranks[p_id].push_back(el_type);

            //====================Hybrid=======================
            for(int k=0;k<nvPerEl;k++)//This works for hexes.
            {
                v_id = Elem2Vert_uniform[el_id][k];
                vertIDs_to_send_to_ranks[p_id].push_back(v_id);
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            
            //====================Hybrid=======================
            not_on_rank++;
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;

            for(int k=0;k<nvPerEl;k++)// looping over the vertices for element "i".
            {
                v_id = Elem2Vert_uniform[el_id][k];
                
                if(m_vertIDs_on_rank.find( v_id ) == m_vertIDs_on_rank.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    m_vertIDs_on_rank.insert(v_id);
                    
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
            
            loc_r_elem.push_back(el_id);
            loc_r_nv_elem.push_back(nvPerEl);
            m_elem2type_on_rank[el_id] = el_type;
            
            on_rank++;
        }
    }

    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);

    std::map<int,std::vector<int> > TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> > part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> > part_tot_recv_elNVs_map;
    
    std::map<int,std::vector<int> >::iterator it;
    
    int n_req_recv;
    int n_req_recv_v;
    int n_req_recv_f;
    int n_req_recv_d;
    std::map<int,std::vector<int> >::iterator totrecv;
    std::map<int,std::vector<int> > part_tot_recv_elTypes_map;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = elms_to_send_to_ranks.begin(); it != elms_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = vertIDs_to_send_to_ranks[it->first].size();

                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, comm);

                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&nvPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*33333+7777, comm);
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&elmtype_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*55555+7777, comm);

                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111, comm, MPI_STATUS_IGNORE);

            std::vector<int>        part_recv_el_id(n_req_recv,0);
            std::vector<int>        part_recv_el_nv(n_req_recv,0);
            std::vector<int>        part_recv_vrt_id(n_req_recv_v,0);
            std::vector<double>     part_recv_vrt_coords(n_req_recv_v*3,0);
            std::vector<int>        part_recv_el_type(n_req_recv,0);

            part_tot_recv_elIDs_map[q]      = part_recv_el_id;

            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_nv[0], n_req_recv, MPI_INT, q, rank*33333+7777, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0],  n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_type[0], n_req_recv, MPI_INT, q, rank*55555+7777, comm, MPI_STATUS_IGNORE);

            part_tot_recv_elTypes_map[q]    = part_recv_el_type;
            TotRecvElement_IDs_v_map[q]     = part_recv_vrt_id;
            part_tot_recv_elIDs_map[q]      = part_recv_el_id;
            part_tot_recv_elNVs_map[q]      = part_recv_el_nv;
        }
    }

    std::vector<int> TotRecvElement_NVs;
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvElement_Types;
    int TotNelem_recv = 0;
    std::vector<int> TotRecvVerts_IDs;
    for(totrecv=part_tot_recv_elIDs_map.begin();totrecv!=part_tot_recv_elIDs_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElement_IDs.push_back(part_tot_recv_elIDs_map[totrecv->first][r]);
            TotRecvElement_NVs.push_back(part_tot_recv_elNVs_map[totrecv->first][r]);
            TotRecvElement_Types.push_back(part_tot_recv_elTypes_map[totrecv->first][r]);
        }

        TotNelem_recv = TotNelem_recv + totrecv->second.size();
    }

    for(totrecv=TotRecvElement_IDs_v_map.begin();totrecv!=TotRecvElement_IDs_v_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvVerts_IDs.push_back(TotRecvElement_IDs_v_map[totrecv->first][r]);
        }
    }

    int Nel_extra   = TotNelem_recv;
    int cnt_v       = 0;
    int cnt_f       = 0;
    int cnt_data    = 0;


    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;

        int nvPerEl = TotRecvElement_NVs[i]; 

        for(int k=0;k<nvPerEl;k++)
        {
            int v_id_n = TotRecvVerts_IDs[cnt_v+k];
            
            if(m_vertIDs_on_rank.find( v_id_n ) == m_vertIDs_on_rank.end()) // add the required unique vertex for current rank.
            {
                m_vertIDs_on_rank.insert(v_id_n);
                
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
        
        cnt_v=cnt_v+nvPerEl;
        on_rank++;
    }


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
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                
                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            //std::cout << rank << " n_reqstd_ids " << n_reqstd_ids << std::endl;
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
        }
    }

    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nverts;
    std::map<int,std::vector<double> > recv_back_verts;
    std::map<int,std::vector<int> > recv_back_verts_ids;
    int n_recv_back;
        
    int nfound=0;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                int dest = it->first;
                //double* vert_send = new double[nv_send*3];
                std::vector<double> vert_send(nv_send*3);
                offset_xcn        = vert_pstate->getOffset(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    if(vertices_i.find(it->second[u])==vertices_i.end())
                    {
                        //std::cout << "NOT FOUND" << std::endl;
                        //std::cout << " it->second[u] " << it->second[u] << " " << dest << " " << rank  << " " << m_Nv_glob << " " << vertices_i.size() << std::endl;
                        nfound++;
                    }
                    else
                    {
                        vert_send[u*3+0]=vertices_i[it->second[u]][0];
                        vert_send[u*3+1]=vertices_i[it->second[u]][1];
                        vert_send[u*3+2]=vertices_i[it->second[u]][2];
                    }
                }
                
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&vert_send.data()[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, comm);
                MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
                //delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<double> recv_back_arr(n_recv_back*3);
            std::vector<int> recv_back_arr_ids(n_recv_back);
            MPI_Recv(&recv_back_arr.data()[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr_ids.data()[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q]     = n_recv_back;
            recv_back_verts[q]      = recv_back_arr;
            recv_back_verts_ids[q]  = recv_back_arr_ids;
        
         }
    }
    
    int vfor = 0;
    std::map<int,std::vector<double> >::iterator it_f;
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


        //std::cout << "giv " << gvid << std::endl;
       
        std::vector<double> V(3,0.0);
        // if(vertices_i.find(gvid)==vertices_i.end())
        // {
        //     std::cout << "Not here " << std::endl;
        // }
        V[0] = vertices_i[gvid][0];
        V[1] = vertices_i[gvid][1];
        V[2] = vertices_i[gvid][2];
        m_LocalVertsMap[gvid]  = V;

        m_LocalV2GlobalV[lvid] = gvid;
        m_GlobalV2LocalV[gvid] = lvid;

        lvid++;
    }
    
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = rank2req_vert[it_f->first][u];
            
            std::vector<double> V(3);
            V[0] = it_f->second[u*3+0];
            V[1] = it_f->second[u*3+1];
            V[2] = it_f->second[u*3+2];
            m_LocalVertsMap[gvid] = V;

            m_LocalV2GlobalV[lvid] = gvid;
            m_GlobalV2LocalV[gvid] = lvid;
           
            m++;
            lvid++;
        }
    }
    
    int nLoc_Verts = m_LocalVertsMap.size();

    int glob_v      = 0;
    int loc_v       = 0;
    int glob_f      = 0;
    int glob_e      = 0;
    double varia_v  = 0.0;
    int ndatas      = 0;

    std::map<int,int> LocElem2Ndata;

    for(m=0;m<loc_r_elem.size();m++)
    {
        
        el_id   = loc_r_elem[m];
        m_ElemSet.insert(el_id);
        nvPerEl = loc_r_nv_elem[m];
        // nfPerEl = loc_r_nf_elem[m];
        // ndatas  = loc_r_Ndata[m];

        // // LocElem2Nv[el_id]    = nvPerEl;
        // // LocElem2Nf[el_id]    = nfPerEl;
        // LocElem2Ndata[el_id] = ndatas;

        // int ndatas = 3;//Elem2Data_uniform[el_id].size();
        // std::vector<double> rowdata(ndatas,0.0);
        // for(int k=0;k<ndatas;k++)
        // {
        //     rowdata[k] = 0.0;//Elem2Data_uniform[el_id][k];
        // }
        
        std::vector<int> tmp_globv;
        std::vector<int> tmp_locv;
        // std::vector<int> tmp_globf;
        // std::vector<int> tmp_globe;
        // std::vector<int> tmp_loce;
        std::vector<double> Pijk(nvPerEl*3);

        for(int p=0;p<nvPerEl;p++)
        {

            glob_v  = Elem2Vert_uniform[el_id][p];
            loc_v   = m_GlobalV2LocalV[glob_v];
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            Pijk[p*3+0]     = m_LocalVertsMap[glob_v][0];
            Pijk[p*3+1]     = m_LocalVertsMap[glob_v][1];
            Pijk[p*3+2]     = m_LocalVertsMap[glob_v][2];

            //m_GlobalVert2Elem[glob_v].push_back(el_id); 
            //globVerts2globElem[glob_v].push_back(el_id);

            //m_Elem2LocalVert[el_id].push_back(loc_v);
            m_Elem2LocalVert[el_id].push_back(loc_v); //globElem2locVerts[el_id].push_back(loc_v);


        }

        std::vector<double> Vijk    = ComputeCentroidCoord(Pijk,nvPerEl);

        //m_Elem2Centroid[el_id]      = Vijk;

        // for(int p=0;p<nfPerEl;p++)
        // {
        //     glob_f = Elem2Face_uniform[el_id][p];
        //     glob_e = Elem2Elem_uniform[el_id][p];

        //     m_globElem2globFaces[el_id].push_back(glob_f);
        //     m_globFace2GlobalElements[glob_f].push_back(el_id);

        //     tmp_globf.push_back(glob_f);
        //     tmp_globe.push_back(glob_e);
        // }

       
        //=======================================
        m_Elem2Vert[el_id]         = tmp_globv;
        // m_Elem2Face[el_id]      = tmp_globf;
        // m_Elem2Elem[el_id]      = tmp_globe;
        // m_Elem2Data[el_id]      = rowdata;

        //=======================================

        tmp_globv.clear();
        tmp_locv.clear();

        /**/
    }
    
    int cnv = 0;
    int cnf = 0;

    for(m=0;m<Nel_extra;m++)
    {
        
        el_id                             = TotRecvElement_IDs[m];
        m_ElemSet.insert(el_id);
        int el_type                       = TotRecvElement_Types[m];
        nvPerEl                           = TotRecvElement_NVs[m];
        m_elem2type_on_rank[el_id]        = el_type;

        std::vector<int> tmp_globv;
        std::vector<int> tmp_locv;

        std::vector<double> Pijk(nvPerEl*3);
        for(int p=0;p<nvPerEl;p++)
        {
            glob_v  = TotRecvVerts_IDs[cnv+p];
            loc_v   = m_GlobalV2LocalV[glob_v];

            m_Elem2LocalVert[el_id].push_back(loc_v);
            // m_GlobalVert2Elem[glob_v].push_back(el_id);

            Pijk[p*3+0] = m_LocalVertsMap[glob_v][0];
            Pijk[p*3+1] = m_LocalVertsMap[glob_v][1];
            Pijk[p*3+2] = m_LocalVertsMap[glob_v][2];

            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            
        }

        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,nvPerEl);

        m_Elem2Centroid[el_id] = Vijk;

        // for(int p=0;p<nfPerEl;p++)
        // {
        //     glob_f = TotRecvFaces_IDs[cnf+p];
        //     glob_e = TotRecvElem_IDs[cnf+p];

        //     m_globElem2globFaces[el_id].push_back(glob_f);
        //     m_globFace2GlobalElements[glob_f].push_back(el_id);

        //     tmp_globf.push_back(glob_f);
        //     tmp_globe.push_back(glob_e);
        // }

        //=======================================
        m_Elem2Vert[el_id]    = tmp_globv;
        // m_Elem2Face[el_id]    = tmp_globf;
        // m_Elem2Elem[el_id]    = tmp_globe;
        // m_Elem2Data[el_id]    = rowdata;
        //=======================================
       

        cnv=cnv+nvPerEl;
        cnf=cnf+nfPerEl;

        tmp_globv.clear();
        tmp_locv.clear();
    }
    

    int nLoc_Elem   = m_Elem2Vert.size();
    vloc            = m_LocalVertsMap.size();

}



void PartObjectLite::ComputeFaceMap(MPI_Comm comm, FaceSetPointer allbcFaces)
{
    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    int fcnt = 0;

    FaceSetPointer m_FaceSetPointer;
    FaceSetPointer m_InteriorFaceSetPointer_tmp;
    
    int ishere = 0;
    int fintid = 0;
    std::vector<std::vector<int> > e2f_map;

    std::map<int,std::vector<int> >::iterator ite;
    int bcf  = 0;
    int intf = 0;
    //std::cout << "m_Elem2Vert " << m_Elem2Vert.size() << m_ElemSet.size() << std::endl;
    for(ite=m_Elem2Vert.begin();ite!=m_Elem2Vert.end();ite++)
    {
        int elid                = ite->first;
        int eltype              = m_elem2type_on_rank[elid];
        std::vector<int> e2n_i  = ite->second;
        int Nv                  = m_Elem2Vert[elid].size();

        std::vector<double> Padj(Nv*3);
        for(int p=0;p<Nv;p++)
        {
            int glob_v   = m_Elem2Vert[elid][p];
            Padj[p*3+0]  = m_LocalVertsMap[glob_v][0];
            Padj[p*3+1]  = m_LocalVertsMap[glob_v][1];
            Padj[p*3+2]  = m_LocalVertsMap[glob_v][2];
        }
        std::vector<double> Vadj  = ComputeCentroidCoord(Padj,Nv);

        switch (eltype) {
        case 10:
            e2f_map = tetra_faces;
            break;
        case 12:
            e2f_map = hex_faces;
            break;
        case 13:
            e2f_map = prism_faces;
            break;
        case 14:    
            e2f_map = pyramid_faces;
            break;
        }
        
        int nfaces = e2f_map.size();
        for(int f=0;f<nfaces;f++)
        {
            fcnt++;
            int nnodes = e2f_map[f].size();
            std::vector<int> face_globvid(nnodes,0);

            for(int n=0;n<nnodes;n++)
            {
                face_globvid[n] = m_Elem2Vert[elid][e2f_map[f][n]];
            }
            
            FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face_globvid));
            // std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            // testInsPointer = m_FaceSetPointer.insert(facePointer);
            
            // This if statement is check to make sure that the face that is being considered is NOT on the boundary.

            FaceSetPointer::iterator FaceFoundInBoundaries = allbcFaces.find(facePointer);

            if(FaceFoundInBoundaries!=allbcFaces.end())
            {
                std::pair<FaceSetPointer::iterator, bool> BoundaryPointerShared = m_OwnedBoundaryFaceSetPointer.insert(facePointer);
                int bFaceID = (*FaceFoundInBoundaries)->GetFaceID();
                (*BoundaryPointerShared.first)->SetFaceID(bFaceID);

                int bFaceRef = (*FaceFoundInBoundaries)->GetFaceRef();
                (*BoundaryPointerShared.first)->SetFaceRef(bFaceRef);

                std::vector<int> fvid                   = (*FaceFoundInBoundaries)->GetEdgeIDs();
                // for(int r=0;r<fvid.size();r++)
                // {
                //     int vid = fvid[r];

                //     if(m_OwnedNonSharedVerts.find(vid)==m_OwnedNonSharedVerts.end())
                //     {
                //         m_OwnedNonSharedVerts.insert(vid);
                //     }
                // }
            }
            else
            {
                std::pair<FaceSetPointer::iterator, bool> testInsPointer2;
                testInsPointer2 = m_InteriorFaceSetPointer_tmp.insert(facePointer);

                if(testInsPointer2.second)
                {
                    (*testInsPointer2.first)->SetFaceID(fintid);
                    (*testInsPointer2.first)->SetFaceLeftElement(elid);
                    (*testInsPointer2.first)->SetFaceRightElement(-1);
                    fintid++;
                }
                else
                {
                    (*testInsPointer2.first)->SetFaceRightElement(elid);
                    intf=intf+1;
                }
            }
        }
    }

    //FaceSetPointer m_FaceSetPointer;
    // update faceid for the boundary faces

    FaceSetPointer::iterator ftit;
    int sharedF = 0;

    std::vector<int> nVrtsPerSharedFace;
    std::vector<int> SharedFaceRankID;   
    std::vector<int> SharedFaceLeftElementID;       
    std::vector<std::vector<int> > sharedFaces;
    int offset_shared_vrts = 0;
    std::vector<int> ownedInteriorFacesNlocs(m_size,0);
    m_ownedInteriorFacesOffsets.resize(m_size);

    for(ftit=m_InteriorFaceSetPointer_tmp.begin();ftit!=m_InteriorFaceSetPointer_tmp.end();ftit++)
    {
        int ofid        = (*ftit)->GetFaceID();
        int rightFaceID = (*ftit)->GetFaceRightElement();
        int leftFaceID  = (*ftit)->GetFaceLeftElement();

        if(rightFaceID==-1) // shared face of partition or interface between boundary layer mesh and tetrahedra
        {
            std::vector<int> fvid                   = (*ftit)->GetEdgeIDs();
            int nvrts_per_face                      = fvid.size();
            SharedFaceRankID.push_back(m_rank);
            SharedFaceLeftElementID.push_back(leftFaceID);
            nVrtsPerSharedFace.push_back(nvrts_per_face);
            sharedFaces.push_back(fvid);
            FaceSharedPtr facePointer               = std::shared_ptr<NekFace>(new NekFace(fvid));
            (*ftit)->SetFaceLeftElement(leftFaceID);

            std::pair<FaceSetPointer::iterator, bool> testInsPointerShared = m_ExternalFaceForRankFaceSetPointer.insert(facePointer);

            offset_shared_vrts                      = offset_shared_vrts+nvrts_per_face;
            sharedF++;
        }
        else
        {
            std::pair<FaceSetPointer::iterator, bool> testInsPointer = m_OwnedInteriorFaceSetPointer.insert((*ftit));
            std::vector<int> fvid                   = (*ftit)->GetEdgeIDs();
            // for(int r=0;r<fvid.size();r++)
            // {
            //     int vid = fvid[r];

            //     if(m_OwnedNonSharedVerts.find(vid)==m_OwnedNonSharedVerts.end())
            //     {
            //         m_OwnedNonSharedVerts.insert(vid);
            //     }
            // }
        }
    }

    // m_ExternalFaceForRankFaceSetPointer hold all the external faces for this particular rank.

    ownedInteriorFacesNlocs[m_rank] = m_OwnedInteriorFaceSetPointer.size();
    m_NownedInteriorFaces.resize(m_size,0);
    MPI_Allreduce(ownedInteriorFacesNlocs.data(),  m_NownedInteriorFaces.data(),   m_size, MPI_INT, MPI_SUM, comm);
    int offsetI = 0;
    for(int i=0;i<m_size;i++)
    {
        m_ownedInteriorFacesOffsets[i] = offsetI;
        offsetI = offsetI + m_NownedInteriorFaces[i];
    }

    int nLocalInterior = m_OwnedInteriorFaceSetPointer.size();
    int nAllInterior = 0;
    MPI_Allreduce(&nLocalInterior,  &nAllInterior,   1, MPI_INT, MPI_SUM, comm);


    //==============================================================================================================================
    //==============================================================================================================================
    //==============================================================================================================================
    // Next step is to figure out all the unique shared faces.
    // This is done by communicating all the shared faces to all the ranks and loop over them to collect all UNIQUE face definitions.
    
    int nSharedFaces                                = nVrtsPerSharedFace.size();
    DistributedParallelState* pstate_sharedFaces    = new DistributedParallelState(nSharedFaces,comm);
    int nAllSharedFaces                             = pstate_sharedFaces->getNel();
    std::vector<int> globalnVrtsPerSharedFace(nAllSharedFaces,0);
    std::vector<int> globalSharedFaceRankID(nAllSharedFaces,0);
    std::vector<int> globalSharedFaceLeftElementID(nAllSharedFaces,0);

    MPI_Allgatherv(&nVrtsPerSharedFace.data()[0],
    nSharedFaces, MPI_INT,
    &globalnVrtsPerSharedFace.data()[0],
    pstate_sharedFaces->getNlocs(),
    pstate_sharedFaces->getOffsets(),
    MPI_INT,comm);

    MPI_Allgatherv(&SharedFaceRankID.data()[0],
    nSharedFaces, MPI_INT,
    &globalSharedFaceRankID.data()[0],
    pstate_sharedFaces->getNlocs(),
    pstate_sharedFaces->getOffsets(),
    MPI_INT,comm);

     MPI_Allgatherv(&SharedFaceLeftElementID.data()[0],
    nSharedFaces, MPI_INT,
    &globalSharedFaceLeftElementID.data()[0],
    pstate_sharedFaces->getNlocs(),
    pstate_sharedFaces->getOffsets(),
    MPI_INT,comm);

    std::vector<int> shared_faces_flatten;
    for (const auto& row : sharedFaces) 
    {
            shared_faces_flatten.insert(shared_faces_flatten.end(), row.begin(), row.end());
    }
    int nSharedFaceVrts                                = shared_faces_flatten.size();
    DistributedParallelState* pstate_sharedFaceVrts    = new DistributedParallelState(nSharedFaceVrts,comm);
    int nVrtsAllSharedFaces                            = pstate_sharedFaceVrts->getNel();

    std::vector<int> globalSharedFaceVrts(nVrtsAllSharedFaces,0);
    MPI_Allgatherv(&shared_faces_flatten.data()[0],
    shared_faces_flatten.size(), MPI_INT,
    &globalSharedFaceVrts.data()[0],
    pstate_sharedFaceVrts->getNlocs(),
    pstate_sharedFaceVrts->getOffsets(),
    MPI_INT,comm);

    int sharedFound = 0;
    int off = 0;
    int itsin = 0;
    int itsnotin = 0;
    int owned_InterFace = 0;
    for(int i=0;i<nAllSharedFaces;i++)
    {
        int nvrts_per_face        = globalnVrtsPerSharedFace[i];
        int rank_id               = globalSharedFaceRankID[i];
        int left_id               = globalSharedFaceLeftElementID[i];

        std::vector<int> face_globvid(nvrts_per_face,0);
        for(int j=0;j<nvrts_per_face;j++)
        {
            face_globvid[j] = globalSharedFaceVrts[off+j];
        }
        off = off + nvrts_per_face;

        FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face_globvid));
        std::pair<FaceSetPointer::iterator, bool> testInsPointer = m_AllSharedFaceSetPointer.insert(facePointer);

        if(testInsPointer.second) // new insertion
        {
            (*testInsPointer.first)->SetFaceRef(rank_id);
            (*testInsPointer.first)->SetFaceLeftRank(rank_id);
            (*testInsPointer.first)->SetFaceLeftElement(left_id);

            // if(rank_id == m_rank)
            // {
            //     owned_InterFace++;
            // }
            // itsnotin++;
        }
        else // face already exists --> these faces are actually shared and do not include the InterFace (Shell) faces.
        {
            // if sharedFace is already in the set, check that its current rank ID is, if rank_id is lower, then update;
            if(rank_id < (*testInsPointer.first)->GetFaceRef())
            {   
                // In this case the SetFaceLeftElement and SetFaceRightElement entries are set to the rank IDS.
                (*testInsPointer.first)->SetFaceRightElement((*testInsPointer.first)->GetFaceLeftElement());
                (*testInsPointer.first)->SetFaceLeftElement(left_id);
                
                (*testInsPointer.first)->SetFaceRightRank((*testInsPointer.first)->GetFaceRef());
                //(*testInsPointer.first)->SetFaceRightRank((*testInsPointer.first)->GetFaceLeftRank());
                (*testInsPointer.first)->SetFaceLeftRank(rank_id);
                
                (*testInsPointer.first)->SetFaceRef(rank_id);
            }
            else
            {
                
                (*testInsPointer.first)->SetFaceLeftElement((*testInsPointer.first)->GetFaceLeftElement());
                (*testInsPointer.first)->SetFaceRightElement(left_id);

                (*testInsPointer.first)->SetFaceLeftRank((*testInsPointer.first)->GetFaceRef());
                (*testInsPointer.first)->SetFaceRightRank(rank_id);
                //std::cout << "get left right " << (*testInsPointer.first)->GetFaceLeftElement() << " " << (*testInsPointer.first)->GetFaceRightElement() << std::endl;  

            }

            std::pair<FaceSetPointer::iterator, bool> ActualSharedFPointer = m_ActualSharedFaceSetPointer.insert((*testInsPointer.first));


            itsin++;
        }
    }

    //==============================================================================================================================
    //==============================================================================================================================
    //==============================================================================================================================
    // std::cout << m_rank << " " << owned_InterFace << std::endl;

    //itsin should be equal to (nAllSharedFaces-nshellFaces)/2
    //std::cout << "nAllSharedFaces " << nAllSharedFaces << std::endl;
    //std::cout << "its in " << itsin << " " << itsnotin << "  " << m_AllSharedFaceSetPointer.size() << " m_ActualSharedFaceSetPointer " << m_ActualSharedFaceSetPointer.size() << std::endl;
    //============================================================================================
    // int shell  = 0;
    // int nshell = 0;
    


    int sharedFOnRank = 0;
    std::vector<int> ownedSharedFacesNlocs(m_size,0);
    m_ownedSharedFacesOffsets.resize(m_size);
    int fon  = 0;
    int nfon = 0;
    int fid  = nAllInterior;
    // int fid  = 0;
    int shared_vrt_id = 0;
    for(ftit=m_ActualSharedFaceSetPointer.begin();ftit!=m_ActualSharedFaceSetPointer.end();ftit++)
    {
        int faceID = fid;
        (*ftit)->SetFaceID(fid);        
        std::vector<int> fvid       = (*ftit)->GetEdgeIDs();
        int fref                    = faceID;
        (*ftit)->SetFaceRef(fref);
        int rankL                   = (*ftit)->GetFaceLeftRank();
        int rankR                   = (*ftit)->GetFaceRightRank();
        int elemL                   = (*ftit)->GetFaceLeftElement();
        int elemR                   = (*ftit)->GetFaceRightElement();
        // std::cout << "This is not set yet  " << fref << std::endl;
        //m_faces4parmmg.push_back(fid);
        
        for(int i=0;i<fvid.size();i++)
        {
            if(m_ActualSharedVerts_g2l.find(fvid[i])==m_ActualSharedVerts_g2l.end())
            {
                m_ActualSharedVerts_g2l[fvid[i]] = shared_vrt_id;
                m_ActualSharedVerts_l2g[shared_vrt_id] = fvid[i];

                shared_vrt_id = shared_vrt_id + 1;
            }
        }

        if(fref == m_rank)
        {
            FaceSharedPtr facePointerNew = std::shared_ptr<NekFace>(new NekFace(fvid));
            facePointerNew->SetFaceRef(fref);
            facePointerNew->SetFaceLeftRank(rankL);
            facePointerNew->SetFaceRightRank(rankR);
            facePointerNew->SetFaceLeftElement(elemL);
            facePointerNew->SetFaceRightElement(elemR);
            facePointerNew->SetFaceID(faceID);
            std::pair<FaceSetPointer::iterator, bool> testInsPointer = m_OwnedSharedFaceSetPointer.insert(facePointerNew);
        }
        fid++;
    }
    
    ownedSharedFacesNlocs[m_rank] = m_OwnedSharedFaceSetPointer.size();
    m_NownedSharedFaces.resize(m_size,0);
    MPI_Allreduce(ownedSharedFacesNlocs.data(),  m_NownedSharedFaces.data(),   m_size, MPI_INT, MPI_SUM, comm);
    int offsetS = 0;
    for(int i=0;i<m_size;i++)
    {
        m_ownedSharedFacesOffsets[i] = offsetS;
        offsetS = offsetS + m_NownedSharedFaces[i];
    }

    int nLocalShared = m_OwnedSharedFaceSetPointer.size();
    int nAllShared = 0;
    MPI_Allreduce(&nLocalShared,  &nAllShared,   1, MPI_INT, MPI_SUM, comm);

    //std::cout << "nAllInterior " << nAllShared+nAllInterior << std::endl;
    // Interior Faces ID    -> Set in this funtion
    // Shared Faces ID      -> Set in this funtion
    // Inter Faces          -> NOT set in this function
    // Boundary Faces       -> NOT set in this function

    //============================================================================================

    int k = 0;
    for(ftit=m_OwnedInteriorFaceSetPointer.begin();ftit!=m_OwnedInteriorFaceSetPointer.end();ftit++)
    {   
        int faceID = m_ownedInteriorFacesOffsets[m_rank] + k;
        
        m_faces4parmmg.push_back(faceID);
        (*ftit)->SetFaceID(faceID);

        std::vector<int> fvid = (*ftit)->GetEdgeIDs();
        for(int r=0;r<fvid.size();r++)
        {
            int vid = fvid[r];

            if(m_ActualSharedVerts_g2l.find(vid)==m_ActualSharedVerts_g2l.end())
            {
                m_OwnedNonSharedVerts.insert(vid);
            }
        }
        // std::vector<int> face = (*ftit)->GetEdgeIDs();
        // FaceSharedPtr FaceNew = std::shared_ptr<NekFace>(new NekFace(face));
        // FaceNew->SetFaceID(faceID);
        // std::pair<FaceSetPointer::iterator, bool> OverallPointer = m_OverallFaceMapOnRank.insert(FaceNew);
        k++;
    }


    for(ftit=m_OwnedBoundaryFaceSetPointer.begin();ftit!=m_OwnedBoundaryFaceSetPointer.end();ftit++)
    {
        std::vector<int> fvid = (*ftit)->GetEdgeIDs();
        for(int r=0;r<fvid.size();r++)
        {
            int vid = fvid[r];

            if(m_ActualSharedVerts_g2l.find(vid)==m_ActualSharedVerts_g2l.end())
            {
                m_OwnedNonSharedVerts.insert(vid);
            }
        }
    }
}

std::map<int,std::vector<double> > PartObjectLite::GatherSharedVertCoordsOnRoot(MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,std::vector<double> > send2root;
    std::map<int,int>::iterator itm;
    for(itm=m_ActualSharedVerts_g2l.begin();itm!=m_ActualSharedVerts_g2l.end();itm++)
    {
        int vid = itm->first;
        
        if(m_LocalVertsMap.find(vid) != m_LocalVertsMap.end())
        {
            send2root[vid] = m_LocalVertsMap[vid];
        }
    }


    std::map<int,std::vector<double> > SharedVertCoordsOnRoot = GatherGlobalMapOnRoot_T(send2root,comm);

    //std::cout << "SharedVertCoordsOnRoot " << SharedVertCoordsOnRoot.size() << " " << m_ActualSharedVerts_g2l.size() << std::endl;

    return SharedVertCoordsOnRoot;

}

std::map<int,std::vector<int> >  PartObjectLite::CommunicateAdjacencyInfoLocalPartition(MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> > Rank2Elem;
    std::map<int,std::vector<int> >::iterator itm;
    int nreq;

    std::set<int> toquery_elids;
    //std::cout << "m_Elem2Elem " << m_Elem2Elem.size() << std::endl;
    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;

        for(int q=0;q<itm->second.size();q++)
        {
            int adjeid = itm->second[q];

            if(m_ElemSet.find(adjeid)==m_ElemSet.end() && adjeid<m_Ne_glob)
            {
                if(toquery_elids.find(adjeid)==toquery_elids.end())
                {
                    toquery_elids.insert(adjeid);
                }
            }
            else
            {
                m_Elem2Rank[adjeid] = rank;
            }
        }
    }

    
    std::vector<int> toquery_elids_vec(toquery_elids.size(),0);
    std::copy(toquery_elids.begin(), 
                toquery_elids.end(), 
                toquery_elids_vec.begin());


    if(rank != 0)
    {
        nreq = toquery_elids_vec.size();
        MPI_Send(&nreq, 1, MPI_INT, 0, rank*100000, comm);
        MPI_Send(&toquery_elids_vec.data()[0], nreq, MPI_INT, 0, rank*1000000, comm);
        
        int nrec_b;
        MPI_Recv(&nrec_b, 1, MPI_INT, 0, rank*250000, comm, MPI_STATUS_IGNORE);
        std::vector<int> pid_2recvback(nrec_b,0);
        std::vector<int> el_2recvback(nrec_b,0);
        MPI_Recv(&el_2recvback[0], nrec_b,  MPI_INT, 0, rank*225000, comm, MPI_STATUS_IGNORE);
        MPI_Recv(&pid_2recvback[0], nrec_b, MPI_INT, 0, rank*20000, comm, MPI_STATUS_IGNORE);

        for(int i=0;i<nrec_b;i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = pid_2recvback[i];

            m_Elem2Rank[el_id] = p_id;

            if(p_id != -1)
            {
                Rank2Elem[p_id].push_back(el_id);
            }
        }
    }
    else if(rank == 0)
    {
        std::vector<int> nrecv_toquery_elids(size-1,0);
        int accul = 0;
        for(int i=1;i<size;i++)
        {
            int dest = i*100000;
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, i, dest, comm, MPI_STATUS_IGNORE);
            nrecv_toquery_elids[i-1] = n_reqstd_ids;
            accul = accul + n_reqstd_ids;
        }

        for(int i=1;i<size;i++)
        {
            int n_reqstd_ids = nrecv_toquery_elids[i-1];
            std::vector<int> QueryOnRoot(n_reqstd_ids,0);
            std::vector<int> pid_2sendback;
            std::vector<int> el_2sendback;
            MPI_Recv(&QueryOnRoot[0], n_reqstd_ids, MPI_INT, i, i*1000000, comm, MPI_STATUS_IGNORE);
            for(int j=0;j<n_reqstd_ids;j++)
            {
                
                if(m_partGlobalRoot.find(QueryOnRoot[j])!=m_partGlobalRoot.end())
                {
                    pid_2sendback.push_back(m_partGlobalRoot[QueryOnRoot[j]]);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
                else
                {                    
                    pid_2sendback.push_back(-1);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
            }
            int send_b = pid_2sendback.size();
            MPI_Send(&send_b, 1, MPI_INT, i, i*250000, comm);
            MPI_Send(&el_2sendback.data()[0], el_2sendback.size(), MPI_INT, i, i*225000, comm);
            MPI_Send(&pid_2sendback.data()[0], pid_2sendback.size(), MPI_INT, i, i*20000, comm);
        }


        for(int i=0;i<toquery_elids_vec.size();i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = -1;
            
            if(m_partGlobalRoot.find(el_id)!=m_partGlobalRoot.end())
            {
                p_id = m_partGlobalRoot[el_id];
                Rank2Elem[p_id].push_back(el_id);
            }
            m_Elem2Rank[el_id] = p_id;
      
        }
    }
    return Rank2Elem;  
}







void PartObjectLite::BuildPMMGCommunicationData(MPI_Comm comm, FaceSetPointer allbcFaces)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    FaceSetPointer::iterator ftit;
    // FaceSetPointer m_ExternalInterFaceFace;
    int f_id = 0;
    int local_fid = 0;
    int f = 0;
    int bcface          = 0;

    for(ftit=m_ExternalFaceForRankFaceSetPointer.begin();ftit!=m_ExternalFaceForRankFaceSetPointer.end();ftit++)
    {
        FaceSetPointer::iterator ActualSharedFPointer = m_ActualSharedFaceSetPointer.find((*ftit));
        (*ftit)->SetFaceID(f_id);
        if(ActualSharedFPointer!=m_ActualSharedFaceSetPointer.end())
        {
            int faceID   = (*ActualSharedFPointer)->GetFaceID();
            int faceref  = (*ActualSharedFPointer)->GetFaceRef();
            int rL      = (*ActualSharedFPointer)->GetFaceLeftRank();
            int rR      = (*ActualSharedFPointer)->GetFaceRightRank();
            // std::cout << "rankie " << rL << " " << rR << std::endl; 
            int EL      = (*ActualSharedFPointer)->GetFaceLeftElement();
            int ER      = (*ActualSharedFPointer)->GetFaceRightElement();
            // m_faces4parmmg.push_back(faceID);
            //std::cout << rL << " L vs R " << rR << " --::-- " << EL << " eL vs eR " << ER <<  std::endl; 
            
            if(m_globShF2locShF.find(faceID)==m_globShF2locShF.end())
            {
                m_faces4parmmg.push_back(faceID);
                m_globShF2locShF[faceID] = f;
                m_locShF2globShF[f] = faceID;
                std::vector<int> facevrts(3,0);
                facevrts[0] = (*ActualSharedFPointer)->GetEdgeIDs()[0];
                facevrts[1] = (*ActualSharedFPointer)->GetEdgeIDs()[1];
                facevrts[2] = (*ActualSharedFPointer)->GetEdgeIDs()[2];

                //================================================================================================
                std::vector<int> face = (*ActualSharedFPointer)->GetEdgeIDs();
                FaceSharedPtr FaceNew = std::shared_ptr<NekFace>(new NekFace(face));
                FaceNew->SetFaceID(faceID);
                FaceNew->SetFaceRef(faceref);
                std::pair<FaceSetPointer::iterator, bool> OverallPointer = m_SharedFacesForRank.insert(FaceNew);
                (*OverallPointer.first)->SetFaceRef(faceref);
                (*OverallPointer.first)->SetFaceID(faceID);
                //================================================================================================
                
                m_sharedFaceVerts[faceID] = facevrts;
                f++;
            }

            if(rL == m_rank)
            {
                m_ColorsFaces[rR].push_back(faceID);
            }
            else if(rR==m_rank)
            {
                m_ColorsFaces[rL].push_back(faceID);
            }
            //externalF2locF[faceID] = local_fid;
            //notinterface++;

            // std::pair<FaceSetPointer::iterator, bool> OverallPointer = m_OverallFaceMapOnRank.insert((*ActualSharedFPointer));

            local_fid++;
        }
        else
        {
            FaceSharedPtr facePointerNew = std::shared_ptr<NekFace>(new NekFace((*ftit)->GetEdgeIDs()));
            std::pair<FaceSetPointer::iterator, bool> testInsPointer = m_ExternalInterFaceFace.insert(facePointerNew);
        }


        // FaceSetPointer::iterator FaceFoundInInterFace = InterFaceFaces.find((*ftit));
        // if(FaceFoundInInterFace!=InterFaceFaces.end())
        // {
        //     interface++;
        //     local_fid++;
        // }

        FaceSetPointer::iterator FaceFoundInBoundaries = allbcFaces.find((*ftit));
        if(FaceFoundInBoundaries!=allbcFaces.end())
        {
            int faceID  = (*ActualSharedFPointer)->GetFaceID();
            // (*FaceFoundInBoundaries)->SetFaceID(bcface);
            m_faces4parmmg.push_back(faceID);
            // bcface++;
        }

        f_id++;
    }

    int nLocalInterior = m_OwnedInteriorFaceSetPointer.size();
    int nAllInterior = 0;
    MPI_Allreduce(&nLocalInterior,  &nAllInterior,   1, MPI_INT, MPI_SUM, comm);
    int FaceOffset = nAllInterior+m_ActualSharedFaceSetPointer.size();

    std::vector<int> external_offset(m_size,0);
    std::vector<int> external_nlocs(m_size,0);
    std::vector<int> external_nlocs_reduce(m_size,0);

    external_nlocs[m_rank] = m_ExternalInterFaceFace.size();
    MPI_Allreduce(external_nlocs.data(),  external_nlocs_reduce.data(),   m_size, MPI_INT, MPI_SUM, comm);
    int offsetInEx = 0;
    for(int i=0;i<m_size;i++)
    {
        external_offset[i] = offsetInEx;
        offsetInEx = offsetInEx + external_nlocs_reduce[i];
    }


    int fid2 = 0;
    for(ftit=m_ExternalInterFaceFace.begin();ftit!=m_ExternalInterFaceFace.end();ftit++)
    {
        int faceID = FaceOffset + external_offset[m_rank] + fid2;
        (*ftit)->SetFaceID(faceID); 
        (*ftit)->SetFaceRef(13); 
         //================================================================================================
        // std::vector<int> face = (*ftit)->GetEdgeIDs();
        // FaceSharedPtr FaceNew = std::shared_ptr<NekFace>(new NekFace(face));
        // FaceNew->SetFaceID(faceID);
        // std::pair<FaceSetPointer::iterator, bool> OverallPointer = m_OverallFaceMapOnRank.insert(FaceNew);
        //================================================================================================
        
        fid2++;
    }

    //std::cout << m_rank << "  "<< m_ExternalInterFaceFace.size() << std::endl;

    // std::cout << "bcface " << bcface << std::endl;

    m_ncomm             = m_ColorsFaces.size();

    m_color_face        = (int *) malloc(m_ncomm*sizeof(int));
    m_ntifc             = (int *) malloc(m_ncomm*sizeof(int));
    
    m_ifc_tria_loc      = (int **)malloc(m_ncomm*sizeof(int *));
    m_ifc_tria_glo      = (int **)malloc(m_ncomm*sizeof(int *));
    
    int icomm           = 0;
    
    std::map<int,std::vector<int> >::iterator itc;
    
    for(itc=m_ColorsFaces.begin();itc!=m_ColorsFaces.end();itc++)
    {
        m_color_face[icomm]     = itc->first;
        m_ntifc[icomm]          =  itc->second.size();
        m_ifc_tria_loc[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));
        m_ifc_tria_glo[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));

        //std::cout << "rank "<< m_rank << " is connected to "  << m_ncomm << " other ranks. #connection = " << icomm << " is rank " << itc->first << " and rank " <<m_rank << " shares " << itc->second.size() << " faces with rank " << itc->first<< std::endl;
        for(int q=0;q<itc->second.size();q++)
        {
            m_ifc_tria_glo[icomm][q] = itc->second[q]+1;
            m_ifc_tria_loc[icomm][q] = m_globShF2locShF[itc->second[q]]+1;
        }

        icomm++;
    }
}





void PartObjectLite::GenerateFace2ElementMap(std::map<int,std::vector<int> > Face2Elem_i,
    MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,std::vector<int> >::iterator itr;

    for(itr=m_Elem2Face.begin();itr!=m_Elem2Face.end();itr++)
    {
        int elid = itr->first;
        int nf   = m_Elem2Face[elid].size();
        int nv   = m_Elem2Vert[elid].size(); 
        int type = m_elem2type_on_rank[elid];

        for(int q = 0; q < m_Elem2Face[elid].size(); q++)
        {
            int fid = m_Elem2Face[elid][q];
            m_Face2Elem[fid].push_back(elid);
        }
    }

    std::map<int,std::vector<int> >::iterator itmiv;
    int tel = 0;

    ParallelState* ife_pstate = new ParallelState(m_Nf_glob,comm);    
    std::map<int,std::vector<int> > rank2req_Faces;
    std::vector<int> new_offsets(size,0);
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }


    for(itmiv=m_Face2Elem.begin();itmiv!=m_Face2Elem.end();itmiv++)
    {
        int faceid   = itmiv->first;
        int numel    = itmiv->second.size();

        if(numel == 1)
        {
            m_Face2Elem[faceid].clear();

            int r = FindRank(new_offsets.data(),size,faceid);

            if(r != rank)
            {
                rank2req_Faces[r].push_back(faceid);
            }
            else
            {
                if(Face2Elem_i.find(faceid)!=Face2Elem_i.end())
                {
                    std::vector<int> newrow(2,0);
                    newrow[0] = Face2Elem_i[faceid][0];
                    newrow[1] = Face2Elem_i[faceid][1];
                    m_Face2Elem[faceid] = newrow;
                }
            }
        }
    }


    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
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
    std::map<int,std::vector<int> > recv_back_face_ids;
    std::map<int,std::vector<int> > recv_back_ife;
    std::map<int,std::vector<int> > recv_back_face_Ne;
    int n_recv_back;

    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send = it->second.size();
                std::vector<int> ife_send;
                std::vector<int> fncol(it->second.size(),0);
                for(int u=0;u<it->second.size();u++)
                {
                    int ncol = 2;
                    ife_send.push_back(Face2Elem_i[it->second[u]][0]);
                    ife_send.push_back(Face2Elem_i[it->second[u]][1]);
                }

                int nfe_send = ife_send.size();

                int dest = it->first;
                MPI_Send(&nfe_send, 1, MPI_INT, dest, 223*6666+1000*dest, comm);
                MPI_Send(&ife_send[0], nfe_send, MPI_INT, dest, 9876*6666+dest*8888, comm);

            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
        {
            int nfe_recv;
            MPI_Recv(&nfe_recv, 1, MPI_INT, q, 223*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_back_ife_arr(nfe_recv);
            MPI_Recv(&recv_back_ife_arr.data()[0], nfe_recv, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_ife[q]        = recv_back_ife_arr;

        }
    }



    std::map<int,std::vector<int> >::iterator iter;
    for(iter=rank2req_Faces.begin();iter!=rank2req_Faces.end();iter++)
    {
        int recvdrank = iter->first;

        int L = iter->second.size();
        int offset = 0;
        int ncol = 0;
        for(int s=0;s<L;s++)
        {
            int face_id  = rank2req_Faces[recvdrank][s];
            ncol     = 2;

            std::vector<int> ife_loc_row(2,0);
            std::vector<int> ifref_loc_row(1,0);
            ife_loc_row[0]   = recv_back_ife[iter->first][offset+0];
            ife_loc_row[1]   = recv_back_ife[iter->first][offset+1];

            m_Face2Elem[face_id] = ife_loc_row;

            offset = offset + ncol;
        }
    }
}



std::map<int,std::vector<double> > PartObjectLite::RedistributeVertexDataForPartition(MPI_Comm comm, int m_Nv_glob, std::map<int, std::vector<double> > t_hessian)
{
    std::map<int,std::vector<double> > owned_hessian;


    int ndatasize = t_hessian.begin()->second.size();
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    ParallelState* vert_pstate = new ParallelState(m_Nv_glob,comm);
    int* new_V_offsets = new int[size];

    for(int i=0;i<size;i++)
    {
        new_V_offsets[i] = vert_pstate->getOffsets()[i]-1;
    }

    std::map<int,std::vector<int> > rank2req_vert;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<double> >::iterator itm;
    for(itm=m_LocalVertsMap.begin();itm!=m_LocalVertsMap.end();itm++)
    {
        int vid = itm->first;
        int rv  = FindRank(new_V_offsets,size,vid);

        if (rv!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
        {
            rank2req_vert[rv].push_back(vid); // add the vertex id that needs to be requested from rank r.
        }
        else
        {
            vertIDs_on_rank.push_back(vid);  // add the vertex to list that is already available on rank.
            owned_hessian[vid] = t_hessian[vid];
        }
    }
    
    ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);

    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    std::map<int,std::vector<int> >::iterator it;

    for(int q=0;q<size;q++)
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
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 6547*2+dest*2, comm);

                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6547+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 6547*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
        }
    }

    int offset_xcn = 0;
    std::map<int,int > recv_back_Nverts;
    std::map<int,std::vector<double> > recv_back_data;

    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                //double* vert_send = new double[nv_send*3];
                int ncol = t_hessian[it->second[0]].size();

                std::vector<double> vert_send(nv_send*ncol);
                offset_xcn        = vert_pstate->getOffset(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    for(int s=0;s<t_hessian[it->second[u]].size();s++)
                    {
                        vert_send[u*ncol+s]=t_hessian[it->second[u]][s];
                    }     
                }
                int dest = it->first;
                MPI_Send(&nv_send, 1, MPI_INT, dest, 6547+1000*dest, comm);
                MPI_Send(&vert_send.data()[0], nv_send*ncol, MPI_DOUBLE, dest, 6547+dest*8888, comm);
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
        {
            int n_recv_back;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 6547+1000*rank, comm, MPI_STATUS_IGNORE);
            std::vector<double> recv_back_arr(n_recv_back*ndatasize);
            MPI_Recv(&recv_back_arr.data()[0], n_recv_back*ndatasize, MPI_DOUBLE, q, 6547+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q]      = n_recv_back;
            recv_back_data[q]        = recv_back_arr;
        }
    }


    int vfor = 0;
    std::map<int,std::vector<int> >::iterator it_f;
    for(it_f=rank2req_vert.begin();it_f!=rank2req_vert.end();it_f++)
    {
        int Nv = it_f->second.size();
        for(int s=0;s<Nv;s++)
        {   
            int vid = it_f->second[s];
            std::vector<double> StateVec(ndatasize,0.0);
            for(int p=0;p<ndatasize;p++)
            {
                StateVec[p]=recv_back_data[it_f->first][s*ndatasize+p];
            }
            owned_hessian[vid] = StateVec;
        }
    }

    return owned_hessian;
}

std::map<int,std::vector<int> > PartObjectLite::getElem2LocalVertMap()
{
    return m_Elem2LocalVert;
}
std::set<int> PartObjectLite::getLocalElemSet()
{
    return m_ElemSet;
}

std::map<int,std::vector<int> > PartObjectLite::getElem2VertMap()
{
    return m_Elem2Vert;
}

std::map<int,std::vector<double> > PartObjectLite::getLocalVertsMap()
{
    return m_LocalVertsMap;
}

FaceSetPointer PartObjectLite::getAllSharedAndInterFaceFaceMap()
{
    return m_AllSharedFaceSetPointer;
}

FaceSetPointer PartObjectLite::getOwnedSharedAndInterFaceFaceMap()
{
    return m_OwnedSharedFaceSetPointer;
}

FaceSetPointer PartObjectLite::getOwnedBoundaryFaceFaceMap()
{
    return m_OwnedBoundaryFaceSetPointer;
}

FaceSetPointer PartObjectLite::getExternalInterFaceFaceMap()
{
    return m_ExternalInterFaceFace;
}

FaceSetPointer PartObjectLite::getActualSharedFaceMap()
{
    return m_ActualSharedFaceSetPointer;
}

std::map<int,int> PartObjectLite::getActualSharedVerts_Global2LocalMap()
{
    return m_ActualSharedVerts_g2l;
}

std::map<int,int> PartObjectLite::getActualSharedVerts_Local2GlobalMap()
{
    return m_ActualSharedVerts_l2g;
}

std::set<int> PartObjectLite::getOwnedSharedVertsMap()
{
    return m_OwnedSharedVerts;
}

std::set<int> PartObjectLite::getOwnedNonSharedVertsMap()
{
    return m_OwnedNonSharedVerts;
}

FaceSetPointer PartObjectLite::getOwnedInteriorFaceFaceMap()
{
    return m_OwnedInteriorFaceSetPointer;
}

FaceSetPointer PartObjectLite::getExternalFacesForRankFaceMap()
{
    return m_ExternalFaceForRankFaceSetPointer;
}

FaceSetPointer PartObjectLite::getSharedFacesForRankMap()
{
    return m_SharedFacesForRank;
}

std::vector<int> PartObjectLite::getownedInteriorFacesOffsets()
{
    return m_ownedInteriorFacesOffsets;
}
std::vector<int> PartObjectLite::getownedInteriorFacesNlocs()
{
    return m_NownedInteriorFaces;
}

std::vector<int> PartObjectLite::getownedSharedFacesOffsets()
{
    return m_ownedSharedFacesOffsets;
}
std::vector<int> PartObjectLite::getownedSharedFacesNlocs()
{
    return m_NownedSharedFaces;
}

std::map<int,int> PartObjectLite::getLocalVert2GlobalVert()
{
    return m_LocalV2GlobalV;
}

std::map<int,int> PartObjectLite::getGlobalVert2LocalVert()
{
    return m_GlobalV2LocalV;
}

int** PartObjectLite::getParMMGCommFace2GlobalVertMap()
{
    return m_ifc_tria_glo;
}

int** PartObjectLite::getParMMGCommFace2LocalVertMap()
{
    return m_ifc_tria_loc;
}

int* PartObjectLite::getParMMGCommColorFace()
{
    return m_color_face;
}

int* PartObjectLite::getParMMGCommNFacesPerColor()
{
    return m_ntifc;
}
int PartObjectLite::getParMMGNComm()
{
    return m_ncomm;
}

std::vector<int> PartObjectLite::getFace4ParMMG()
{
    return m_faces4parmmg;
}

std::map<int,int> PartObjectLite::getLocalSharedFace2GlobalSharedFace()
{
    return m_locShF2globShF;
}
std::map<int,int> PartObjectLite::getGlobalSharedFace2LocalSharedFace()
{   
    return m_globShF2locShF;
}

std::map<int,std::vector<int> > PartObjectLite::getSharedFaceMap()
{
    return m_sharedFaceVerts;
}

PartObjectLite::~PartObjectLite()
{

}


/*
std::map<int, std::vector<int> > PartObjectLite::getAdjacentElementLayer(std::map<int,std::vector<double> > xcn, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);
    ParallelState* xcn_pstate = new ParallelState(Nv_glob,comm);
    int* new_V_offsets = new int[size];
    for(int i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }
    std::map<int, std::vector<int> > adjEl2Face;

    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int e_id;
    int r;
        int itel = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    int ff   = 0;
    // std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_vert;
    // std::map<int,std::vector<int> > rank2req_face;
    std::map<int, std::vector<int> > adj_el;

    std::map<int,std::vector<int> > Rank2Elem = CommunicateAdjacencyInfoExtendedPartition(comm);
    std::map<int,std::vector<int> >::iterator iter;


    for(iter = Rank2Elem.begin();iter != Rank2Elem.end(); iter++)
    {
        int ra = iter->first;
        
        if(m_Rank2ReqElem.find(ra)==m_Rank2ReqElem.end())
        {
            std::set<int> el_list;
            for(int q=0;q<iter->second.size();q++)
            {
                el_list.insert(iter->second[q]);
            }
            m_Rank2ReqElem[ra] = el_list;
        }
        if(m_Rank2ReqElem.find(ra)!=m_Rank2ReqElem.end())
        {
            for(int q=0;q<iter->second.size();q++)
            {
                m_Rank2ReqElem[ra].insert(iter->second[q]);
            }
        }
    }

    //std::cout << "Funished " << adj_elements.size() << std::endl;
    // adj_elements = adjid2rank;
    // std::cout << "rank " << rank << " adjid2rank " << adjid2rank.size() << std::endl;

    // if(rank == 1)
    // {
    //     std::map<int,int>::iterator itc;
    //     std::cout << "Should ivert right " << std::endl;
    //     for(itc=adjid2rank.begin();itc!=adjid2rank.end();itc++)
    //     {
    //         std::cout << itc->first << " " << itc->second << std::endl;
    //     }
    // }

    ScheduleObj* adj_schedule = DoScheduling(Rank2Elem, comm);
    std::map<int,std::vector<int> > reqstd_adj_ids_per_rank;
    std::map<int,std::vector<int> >::iterator it;
    int n_reqstd_adj_ids;


    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = Rank2Elem.begin(); it != Rank2Elem.end(); it++)
            {
                int n_req_adj_el           = it->second.size();
                int dest                   = it->first;
                //std::cout << dest << " " << adj_elements.size() << std::endl;
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req_adj_el, 1, MPI_INT, dest, 9876000+10*dest, comm);
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req_adj_el, MPI_INT, dest, 9876000*2+dest*2, comm);
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


    //std::cout << reqstd_adj_ids_per_rank.size() << " " << rank << std::endl;

   
    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,std::vector<int> > send_adj_verts_IDs;
    // std::map<int,std::vector<int> > send_adj_faces_IDs;
    std::map<int,std::vector<int> > send_adj_element_IDs;

    std::map<int,std::vector<int> > send_adj_NvertsPel;
    // std::map<int,std::vector<int> > send_adj_NfacesPel;
    std::map<int,std::vector<int> > send_adj_NdatasPel;
    std::map<int,std::vector<double> > send_adj_data;
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

            int nvPerEl = m_Elem2Vert[adj_id].size();
            // int nfPerEl = m_Elem2Face[adj_id].size();
            int ndata   = m_Elem2Data[adj_id].size();
            //Array<double>* Var  = LocAndAdjElemVaria[adj_id];
            
            send_adj_NvertsPel[dest].push_back(nvPerEl);
            // send_adj_NfacesPel[dest].push_back(nfPerEl);
            send_adj_NdatasPel[dest].push_back(ndata);

            for(int k=0;k<ndata;k++)
            {
                double dataV = m_Elem2Data[adj_id][k];
                send_adj_data[dest].push_back(dataV);
            }

            //send_adj_VarPel[dest].push_back();
            
            for(int k=0;k<nvPerEl;k++)
            {
                v_id = m_Elem2Vert[adj_id][k];
                send_adj_verts_IDs[dest].push_back(v_id);
            }
            // for(int k=0;k<nfPerEl;k++)
            // {
            //     f_id = m_Elem2Face[adj_id][k];
            //     send_adj_faces_IDs[dest].push_back(f_id);
            //     e_id = m_Elem2Elem[adj_id][k];
            //     send_adj_element_IDs[dest].push_back(e_id);
            // }
        }

        TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn   = 0;
    std::map<int,int  > recv_adj_back_Nverts;
    std::map<int,std::vector<int> > recv_adj_back_verts_ids;
    std::map<int,int  > recv_adj_back_Nfaces;
    // std::map<int,std::vector<int> > recv_adj_back_faces_ids;
    // std::map<int,std::vector<int> > recv_adj_back_element_ids;

    //std::map<int,int > recv_adj_back_Nrhos;
    std::map<int,std::vector<int>  > recv_adj_NvPel;
    // std::map<int,std::vector<int>  > recv_adj_NfPel;
    std::map<int,std::vector<int>  > recv_adj_NdataPel;
    std::map<int,std::vector<double>  >recv_adj_NVarPel;
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
                MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 19999*9876+dest*8888,comm);
                
                int NnvPel    = send_adj_NvertsPel[it->first].size();
                // int NnfPel    = send_adj_NfacesPel[it->first].size();
                int NndataPel = send_adj_NdatasPel[it->first].size();
                int NVarPel   = send_adj_data[it->first].size();
                MPI_Send(&NnvPel,  1, MPI_INT, dest, 98764444+5000*dest, comm);
                // MPI_Send(&NnfPel,  1, MPI_INT, dest, 98764444-5000*dest, comm);
                MPI_Send(&NndataPel,  1, MPI_INT, dest, 98765555-5000*dest, comm);

                MPI_Send(&NVarPel, 1, MPI_INT, dest, 98766666-5000*dest, comm);

                MPI_Send(&send_adj_NvertsPel[it->first][0], NnvPel, MPI_INT, dest, 98364444+15000*dest, comm);
                // MPI_Send(&send_adj_NfacesPel[it->first][0], NnfPel, MPI_INT, dest, 98364444-15000*dest, comm);
                MPI_Send(&send_adj_NdatasPel[it->first][0], NndataPel, MPI_INT, dest, 98364444-15000*dest, comm);

                MPI_Send(&send_adj_data[it->first][0], NVarPel, MPI_DOUBLE, dest, 98366666-15000*dest, comm);

                // int nf_adj_send = send_adj_faces_IDs[it->first].size();
                // MPI_Send(&nf_adj_send, 1, MPI_INT, dest, 3333*9876+dest*8888,comm);
                // MPI_Send(&send_adj_faces_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 2222*9876+dest*8888,comm);
                // MPI_Send(&send_adj_element_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 5555*9876+dest*8888,comm);

                //int n_adj_rhos = send_adj_rhos[it->first].size();
                //MPI_Send(&n_adj_rhos, 1, MPI_INT, dest, 4444*9876+dest*8888,comm);
                //MPI_Send(&send_adj_rhos[it->first][0], n_adj_rhos, MPI_DOUBLE, dest, 5555*9876+dest*8888,comm);


            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 98760000+1000*rank, comm, MPI_STATUS_IGNORE);
            //int* recv_adj_back_arr_ids = new int[n_adj_vert_recv_back];
            std::vector<int> recv_adj_back_arr_ids(n_adj_vert_recv_back);
            MPI_Recv(&recv_adj_back_arr_ids.data()[0], n_adj_vert_recv_back, MPI_INT, q, 19999*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            
            int NnvPel_recv_back,NnfPel_recv_back,NndataPel_recv_back,NVarPel_recv_back;
            MPI_Recv(&NnvPel_recv_back, 1, MPI_INT, q, 98764444+5000*rank, comm, MPI_STATUS_IGNORE);
            // MPI_Recv(&NnfPel_recv_back, 1, MPI_INT, q, 98764444-5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NndataPel_recv_back, 1, MPI_INT, q, 98765555-5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NVarPel_recv_back, 1, MPI_INT, q, 98766666-5000*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> Nnv_RB(NnvPel_recv_back);
            // std::vector<int> Nnf_RB(NnfPel_recv_back);
            std::vector<int> Nndata_RB(NndataPel_recv_back);
            std::vector<double> NnVar_RB(NVarPel_recv_back);

            MPI_Recv(&Nnv_RB[0], NnvPel_recv_back, MPI_INT, q, 98364444+15000*rank, comm, MPI_STATUS_IGNORE);
            // MPI_Recv(&Nnf_RB[0], NnfPel_recv_back, MPI_INT, q, 98364444-15000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&Nndata_RB[0], NndataPel_recv_back, MPI_INT, q, 98364444-15000*rank, comm, MPI_STATUS_IGNORE);

            MPI_Recv(&NnVar_RB[0], NVarPel_recv_back, MPI_DOUBLE, q, 98366666-15000*rank, comm, MPI_STATUS_IGNORE);

            // int n_adj_face_recv_back;
            // MPI_Recv(&n_adj_face_recv_back, 1, MPI_INT, q, 3333*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            // //int* recv_adj_back_arr_face_ids = new int[n_adj_face_recv_back];
            // std::vector<int> recv_adj_back_arr_face_ids(n_adj_face_recv_back);
            // std::vector<int> recv_adj_back_arr_element_ids(n_adj_face_recv_back);

            // MPI_Recv(&recv_adj_back_arr_face_ids.data()[0], n_adj_face_recv_back, MPI_INT, q, 2222*9876+rank*8888, comm,   MPI_STATUS_IGNORE);
            // MPI_Recv(&recv_adj_back_arr_element_ids.data()[0], n_adj_face_recv_back, MPI_INT, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);

            //int n_adj_rho_recv_back;
            //MPI_Recv(&n_adj_rho_recv_back, 1, MPI_INT, q, 4444*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            //double* recv_adj_back_arr_rho = new double[n_adj_rho_recv_back];
            //MPI_Recv(&recv_adj_back_arr_rho[0], n_adj_rho_recv_back, MPI_DOUBLE, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);


            recv_adj_back_Nverts[q]     = n_adj_vert_recv_back;
            recv_adj_back_verts_ids[q]  = recv_adj_back_arr_ids;

            // recv_adj_back_Nfaces[q]     = n_adj_face_recv_back;
            // recv_adj_back_faces_ids[q]  = recv_adj_back_arr_face_ids;
            // recv_adj_back_element_ids[q]  = recv_adj_back_arr_element_ids;

            recv_adj_NvPel[q] = Nnv_RB;
            // recv_adj_NfPel[q] = Nnf_RB;
            recv_adj_NdataPel[q] = Nndata_RB;
            recv_adj_NVarPel[q] = NnVar_RB;
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

    // std::vector<int> adj_faces;
    // std::vector<int> adj_element;
    // for(itm=recv_adj_back_Nfaces.begin();itm!=recv_adj_back_Nfaces.end();itm++)
    // {
    //     TotNface_adj_recv = TotNface_adj_recv+itm->second;
    //     for(int i=0;i<itm->second;i++)
    //     {
    //         adj_faces.push_back(recv_adj_back_faces_ids[itm->first][i]);
    //         adj_element.push_back(recv_adj_back_element_ids[itm->first][i]);
    //     }
    // }

    std::vector<int> Rank2Elem_vec;
    std::vector<int> NvPEl_rb;
    //std::vector<int> NfPEl_rb;
    std::vector<int> NdataPEl_rb;
    std::vector<double> NVarPEl_rb;
    std::map<int,std::vector<int> >::iterator itm_el;
    int offvvv = 0;
   
    for(itm_el=Rank2Elem.begin();itm_el!=Rank2Elem.end();itm_el++)
    {
        
        for(int i=0;i<itm_el->second.size();i++)
        {
            Rank2Elem_vec.push_back(Rank2Elem[itm_el->first][i]);
            int Nv    = recv_adj_NvPel[itm_el->first][i];
            //int Nf    = recv_adj_NfPel[itm_el->first][i];
            int ndata = recv_adj_NdataPel[itm_el->first][i];
            for(int q=0;q<ndata;q++)
            {
                NVarPEl_rb.push_back(recv_adj_NVarPel[itm_el->first][i*ndata+q]);
            }

            NvPEl_rb.push_back(Nv);
            // NfPEl_rb.push_back(Nf);
            NdataPEl_rb.push_back(ndata);
            
            offvvv=offvvv+Nv;

        }
    }
    
    for(int i=0;i<adj_verts.size();i++)
    {
        int v_id_n = adj_verts[i];
        r = FindRank(new_V_offsets,size,v_id_n);

        if(m_vertIDs_on_rank.find( v_id_n ) == m_vertIDs_on_rank.end()) // add the required unique vertex for current rank.
        {
            m_vertIDs_on_rank.insert(v_id_n);
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
    
    // for(int i=0;i<adj_faces.size();i++)
    // {
    //     int f_id_n = adj_faces[i];

    //     if(m_faceIDs_on_rank.find( f_id_n ) == m_faceIDs_on_rank.end()) // add the required unique vertex for current rank.
    //     {
    //         m_faceIDs_on_rank.insert(f_id_n);
    //         faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
    //         floc_tmp++;

    //         if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
    //         {
    //             rank2req_face[r].push_back(f_id_n); // add vertex to rank2req_vert map.
    //         }
    //         else
    //         {
    //             faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
    //             floc_tmp++;
    //         }
    //     }
        
    // }
    
    

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
                   MPI_Send(&n_req, 1, MPI_INT, dest, 6547+10*dest, comm);
                   //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                   MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 6547*2+dest*2, comm);

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
       std::map<int,std::vector<double> > recv_back_verts;
       std::map<int,std::vector<int> > recv_back_verts_ids;
       int n_recv_back;
    
    
    

        
    
       for(q=0;q<size;q++)
       {
           if(rank == q)
           {
               for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
               {
                   int nv_send = it->second.size();
                   //double* vert_send = new double[nv_send*3];
                   std::vector<double> vert_send(nv_send*3);
                   offset_xcn        = xcn_pstate->getOffset(rank);
                   for(int u=0;u<it->second.size();u++)
                   {
                       // vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                       // vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                       // vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);

                        vert_send[u*3+0]=xcn[it->second[u]][0];
                        vert_send[u*3+1]=xcn[it->second[u]][1];
                        vert_send[u*3+2]=xcn[it->second[u]][2];
                   }

                   int dest = it->first;
                   MPI_Send(&nv_send, 1, MPI_INT, dest, 6547+1000*dest, comm);
                   // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);

                   MPI_Send(&vert_send.data()[0], nv_send*3, MPI_DOUBLE, dest, 6547+dest*8888, comm);
                   MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 8888*6547+dest*8888,comm);

                   //delete[] vert_send;
               }
           }
           if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
            {
               MPI_Recv(&n_recv_back, 1, MPI_INT, q, 6547+1000*rank, comm, MPI_STATUS_IGNORE);

               //double* recv_back_arr = new double[n_recv_back*3];
               std::vector<double> recv_back_arr(n_recv_back*3);
               //int* recv_back_arr_ids = new int[n_recv_back];
               std::vector<int> recv_back_arr_ids(n_recv_back);
               //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr.data()[0], n_recv_back*3, MPI_DOUBLE, q, 6547+rank*8888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr_ids.data()[0], n_recv_back, MPI_INT, q, 8888*6547+rank*8888, comm, MPI_STATUS_IGNORE);

               recv_back_Nverts[q]     = n_recv_back;
               recv_back_verts[q]      = recv_back_arr;
               recv_back_verts_ids[q]  = recv_back_arr_ids;

               }
       }


       
       int vfor = 0;
       std::map<int,std::vector<double> >::iterator it_f;
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

           if(m_GlovalV2LocalV.find(gvid)==m_GlovalV2LocalV.end())
           {
               
               std::vector<double> V(3);

               V[0] = xcn[gvid][0];
               V[1] = xcn[gvid][1];
               V[2] = xcn[gvid][2];
               
               m_LocalVertsMap[gvid] = V;

               lvid++;
           }
       }

       //int o = 3*vloc_tmp;
       int u = 0;
       for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
       {
           int Nv = recv_back_Nverts[it_f->first];

           for(u=0;u<Nv;u++)
           {
               gvid = recv_back_verts_ids[it_f->first][u];
               
               if(m_GlovalV2LocalV.find(gvid)==m_GlovalV2LocalV.end())
               {
                   
                   std::vector<double> V(3);
                   V[0] = it_f->second[u*3+0];
                   V[1] = it_f->second[u*3+1];
                   V[2] = it_f->second[u*3+2];

                   //LocalVerts.push_back(V);
                   m_LocalVertsMap[gvid] = V;
                //    o_lvertex2gvertex[lvid]=gvid;
                //    m_GlovalV2LocalV[gvid]=lvid;
                   
                   lvid++;
               }
            }
       }

    int nLoc_Verts = m_LocalVertsMap.size();

    // ================================== Faces on Rank =========================================


    int cnv = 0;
    int cnf = 0;
    int idsave = 0;
    //double rho_v;
    int loc_v;
    int glob_f;
    int loc_f;
    int glob_v;
    int glob_e=0;
    int before = m_Elem2Elem.size();

    double varia_v = 0.0;
    std::set<int> LocAdjElemSet;
    std::vector<int> adjElLayer(3*itel);
    int offvv = 0;
    for(int m=0;m<Rank2Elem_vec.size();m++)
    {
        el_id         = Rank2Elem_vec[m];
        
        int Nv        = NvPEl_rb[m];
        // int Nf        = NfPEl_rb[m];
        int ndatas    = NdataPEl_rb[m];

        std::vector<int> tmp_globv(Nv,0);
        std::vector<int> tmp_locv(Nv,0);
        // std::vector<int> tmp_globf(Nf,0);
        // std::vector<int> tmp_globe(Nf,0);

        std::vector<double> datarow(ndatas,0.0);
        for(int q=0;q<ndatas;q++)
        {
            datarow[q] = NVarPEl_rb[m*ndatas+q];
        }
    
        
        std::vector<double> Padj(Nv*3);
        for(int p=0;p<Nv;p++)
        {
            glob_v = adj_verts[offvv+p];
            loc_v  = m_GlovalV2LocalV[glob_v];

            Padj[p*3+0] = m_LocalVertsMap[glob_v][0];
            Padj[p*3+1] = m_LocalVertsMap[glob_v][1];
            Padj[p*3+2] = m_LocalVertsMap[glob_v][2];

            tmp_globv[p] = glob_v;
            tmp_locv[p] = loc_v;

            m_Elem2LocalVert[el_id].push_back(loc_v);//globElem2locVerts[el_id].push_back(loc_v);
            cnv++;
            
        }
        std::vector<double> Vadj  = ComputeCentroidCoord(Padj,Nv);
        m_Elem2Centroid[el_id] = Vadj;

        // for(int p=0;p<Nf;p++)
        // {
        //     glob_f = adj_faces[cnf];
        //     glob_e = adj_element[cnf];
        //     tmp_globf[p] = glob_f;
        //     tmp_globe[p] = glob_e;
        //     m_globFace2GlobalElements[glob_f].push_back(el_id);
        //     cnf++;
        // }
    
        //=======================================
        m_Elem2Vert[el_id]      = tmp_globv;
        m_Elem2Face[el_id]      = tmp_globf;
        m_Elem2Elem[el_id]      = tmp_globe;
        m_Elem2Data[el_id]      = datarow;
        //=======================================
        adjEl2Face[el_id]       = tmp_globf;

        offvv = offvv+Nv;
        tmp_locv.clear();
    }

    delete[] new_V_offsets;

    reqstd_ids_per_rank.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    
    NvPEl_rb.clear();
    NfPEl_rb.clear();
    
    rank2req_vert.clear();


    return adjEl2Face;

}*/