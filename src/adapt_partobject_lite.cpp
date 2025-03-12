


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
    /**/
    std::cout << "m_Elem2Vert  " << m_Elem2Vert.size() << " on rank " << rank << " " << nfound << std::endl;


}