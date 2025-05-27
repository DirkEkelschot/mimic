


#include "adapt_partobject_lite.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"


PartObjectLite::PartObjectLite(std::map<int,std::vector<int> > Elem2Vert_uniform,
                               std::map<int,std::vector<double> > vertices_i,
                               std::map<int,int> eltype_map,
                               std::vector<int> elTypes,
                               int nE,
                               int nV,
                               MPI_Comm comm)
                               
{

    m_Ne_glob = nE;
    m_Nv_glob = nV;

    DeterminePartitionLayout(Elem2Vert_uniform,elTypes,comm);

    DetermineElement2ProcMap(Elem2Vert_uniform, 
                             vertices_i,
                             eltype_map, 
                             comm);


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
            // loc_v   = m_GlobalV2LocalV[glob_v];

            // if(glob_v>vert_pstate->getNel())
            // {
            //     std::cout << "Nel On error " << glob_v << std::endl;
            // }
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            Pijk[p*3+0]     = m_LocalVertsMap[glob_v][0];
            Pijk[p*3+1]     = m_LocalVertsMap[glob_v][1];
            Pijk[p*3+2]     = m_LocalVertsMap[glob_v][2];


            m_GlobalVert2Elem[glob_v].push_back(el_id); //globVerts2globElem[glob_v].push_back(el_id);
            // m_Elem2LocalVert[el_id].push_back(loc_v); //globElem2locVerts[el_id].push_back(loc_v);
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
        m_Elem2Vert[el_id]      = tmp_globv;
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

            // m_Elem2LocalVert[el_id].push_back(loc_v);
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
    // floc            = cnf;
    
    std::cout << "m_Elem2Vert  " << m_Elem2Vert.size() << " on rank " << rank << " " << nfound << std::endl;
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