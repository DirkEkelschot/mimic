#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted
void merge_indexed(int data[], const int offsets[], size_t index_begin, size_t index_end)
{
    if (index_end - index_begin > 1) {
        auto index_middle = index_begin + (index_end - index_begin) / 2;
        merge_indexed(data, offsets, index_begin, index_middle);
        merge_indexed(data, offsets, index_middle, index_end);
        std::inplace_merge(&data[offsets[index_begin]], &data[offsets[index_middle]], &data[offsets[index_end]]);
    }
}

struct PartMesh
{
    std::map<int,int> LocalVert2GlobalVert;
    std::map<int,int> GlobalVert2LocalVert;
    std::vector<Vert*> LocalVerts;
    std::map<int, std::vector<int> > ien_part;
    std::map<int, std::vector<int> > ien_part_glob;
    std::map<int, std::vector<int> > ief_part;
};


struct TetrahedraMesh
{
    Array<int>* ien_part_tetra;
    Array<int>* ief_part_tetra;
    
    Array<int>* ien_part_hybrid;
    
};


PartMesh* CommunicateTetrahedra(int nglob, Array<int>* part,
                           Array<int>* new_ien,
                           Array<int>* new_ien_or,
                           Array<int>* new_ief,
                           ParArray<double>* xcn,
                           ParallelState* xcn_pstate,
                           ParArray<int>* ife,
                           ParallelState* ife_pstate,
                           MPI_Comm comm)
{
    
    PartMesh* pm = new PartMesh;
    ParallelState* ien_pstate               = new ParallelState(nglob,comm);

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
    int v_id_o;
    Vert V;
    std::vector<Vert> part_verts;
    std::vector<std::vector<int> >  part_elem2verts;
    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > OriVertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_vertOri;
    std::map<int,std::vector<int> > rank2req_face;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::vector<int> vertIDs_on_rank_Ori;
    std::vector<int> part_v;
    std::vector<int> loc_r_elem;
    
    std::set<int> unique_vertIDs_on_rank_set;
    std::set<int> unique_faceIDs_on_rank_set;
    
    int r     = 0;
    int lv_id = 0;
    int lf_id = 0;
    int f_id  = 0;

    int xcn_o = xcn->getOffset(rank);
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
        el_id   = ien_pstate->getOffsets()[rank]+i;      
        
        nvPerEl = 4;
        nfPerEl = 4;
        
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.

            //====================Hybrid=======================
            for(int k=0;k<4;k++)
            {
                v_id   = new_ien->getVal(i,k);
                v_id_o = new_ien_or->getVal(i,k);
                vertIDs_to_send_to_ranks[p_id].push_back(v_id);
                OriVertIDs_to_send_to_ranks[p_id].push_back(v_id_o);
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            
            
            
            for(int k=0;k<4;k++)
            {
                f_id = new_ief->getVal(i,k);
                faceIDs_to_send_to_ranks[p_id].push_back(f_id);
            }
            //====================Hybrid=======================
            not_on_rank++;
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;
            std::vector<int> nodes(4);
            std::vector<int> nodes_glob(4);
            std::vector<int> faces(4);

            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                v_id   = new_ien->getVal(i,k);
                v_id_o = new_ien_or->getVal(i,k);
                nodes[k] = v_id;
                nodes_glob[k] = v_id_o;
                
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_V_offsets,size,v_id_o);
                    
                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                        rank2req_vertOri[r].push_back(v_id_o);
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vertIDs_on_rank_Ori.push_back(v_id_o);
//                        if(rank == 2)
//                        {
//                            std::cout << "inside " << v_id_o << " " << r << " " << new_V_offsets[rank] << std::endl;
//                        }
                        
                        vloc_tmp++;
                    }
                    lv_id++;
                }
                
                
                
            }
            
            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                f_id = new_ief->getVal(i,k);
                faces[k] = f_id;
                
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
            //std::cout << "el ID = " << el_id << " wr " << rank << std::endl;
            pm->ien_part[el_id]      = nodes;
            pm->ien_part_glob[el_id] = nodes_glob;
            pm->ief_part[el_id]      = faces;
            loc_r_elem.push_back(el_id);
            on_rank++;
        }
        
        
    }
    
    std::cout << " RANK - " << rank << " " <<vertIDs_on_rank.size() << std::endl;
    
    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);
    std::map<int,std::vector<int> >  part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_ov_map;
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
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&OriVertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 339000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222, comm, MPI_STATUS_IGNORE);

            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_vrt_id(n_req_recv_v);
            std::vector<int>    part_recv_orivrt_id(n_req_recv_v);
            std::vector<int>    part_recv_face_id(n_req_recv_f);
            
            MPI_Recv(&part_recv_el_id[0],   n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0],  n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orivrt_id[0],  n_req_recv_v, MPI_INT, q, 339000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            TotRecvElement_IDs_ov_map[q] = part_recv_orivrt_id;
            TotRecvElement_IDs_v_map[q] = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q] = part_recv_face_id;
            part_tot_recv_elIDs_map[q]  = part_recv_el_id;
        }
    }
    
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvVerts_IDs;
    std::vector<int> TotRecvOriVerts_IDs;
    std::vector<int> TotRecvFaces_IDs;

    std::map<int,std::vector<int> >::iterator totrecv;
    //unpack the element IDs and their corresponding variable values.
    int TotNelem_recv = 0;
    for(totrecv=part_tot_recv_elIDs_map.begin();totrecv!=part_tot_recv_elIDs_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElement_IDs.push_back(part_tot_recv_elIDs_map[totrecv->first][r]);
        }
        TotNelem_recv = TotNelem_recv + totrecv->second.size();
    }
    //unpack the vertex IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_v_map.begin();totrecv!=TotRecvElement_IDs_v_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvVerts_IDs.push_back(TotRecvElement_IDs_v_map[totrecv->first][r]);
            TotRecvOriVerts_IDs.push_back(TotRecvElement_IDs_ov_map[totrecv->first][r]);
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
        std::vector<int> nodes(4);
        std::vector<int> nodes_glob(4);
        std::vector<int> faces(4);
        
        int elID = TotRecvElement_IDs[i];
        for(int k=0;k<4;k++)
        {
            int v_id_n = TotRecvVerts_IDs[cnt_v+k];
            int v_id_o_n = TotRecvOriVerts_IDs[cnt_v+k];
            
            nodes[k] = v_id_n;
            nodes_glob[k] = v_id_o_n;
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_V_offsets, size, v_id_o_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_vert[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                    rank2req_vertOri[r].push_back(v_id_o_n);
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    vertIDs_on_rank_Ori.push_back(v_id_o_n);
                    vloc_tmp++;
                }
                lv_id++;
            }
        }
        for(int k=0;k<4;k++)// looping over the vertices for element "i".
        {
            int f_id_n = TotRecvVerts_IDs[cnt_v+k];
            faces[k] = f_id_n;
            
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
        
        cnt_v=cnt_v+4;
        cnt_f=cnt_f+4;
        
        pm->ien_part[elID] = nodes;
        pm->ien_part_glob[elID] = nodes_glob;
        pm->ief_part[elID] = faces;
        
        on_rank++;
        
    }
    
    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    
    
    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    int m = 0;
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Ori_ids_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_vert.begin(); it != rank2req_vert.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                
                //	MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //	MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                MPI_Send(&rank2req_vertOri[it->first][0], n_req, MPI_INT, dest, 2229876*2+dest*2, comm);

                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_Ori_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_reqstd_Ori_ids[0], n_reqstd_ids, MPI_INT, q, 2229876*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
            reqstd_Ori_ids_per_rank[q] = recv_reqstd_Ori_ids;
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
            for (it = reqstd_Ori_ids_per_rank.begin(); it != reqstd_Ori_ids_per_rank.end(); it++)
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
                MPI_Send(&reqstd_ids_per_rank[it->first][0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
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
    int gvid_gl = 0;

    for(m=0;m<vloc_tmp;m++)
    {
        gvid = vertIDs_on_rank[m];
        gvid_gl = vertIDs_on_rank_Ori[m];
        
        Vert* V = new Vert;
        
        V->x = xcn->getVal(gvid_gl-xcn_o,0);
        V->y = xcn->getVal(gvid_gl-xcn_o,1);
        V->z = xcn->getVal(gvid_gl-xcn_o,2);
            
        pm->LocalVerts.push_back(V);
        pm->LocalVert2GlobalVert[lvid] = gvid;
        pm->GlobalVert2LocalVert[gvid] = lvid;
        lvid++;
    }
    
    m = 0;
    
    std::cout << "# verts intermediate " << rank << " " << pm->LocalVerts.size() << " vloc " << vloc_tmp << std::endl;
    
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
            
            pm->LocalVerts.push_back(V);
            
            pm->LocalVert2GlobalVert[lvid]=gvid;
            pm->GlobalVert2LocalVert[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }


    int nLoc_Verts = pm->LocalVerts.size();
    
    
    std::cout << "# verts final " << rank << " " << pm->LocalVerts.size() << " vloc " << vloc_tmp << std::endl;
    // ================================== Faces on Rank =========================================
    
//    int lfid = 0;
//    int gfid = 0;
//    for(m=0;m<floc_tmp;m++)
//    {
//        gfid = faceIDs_on_rank[m];
//
//        LocalFace2GlobalFace[lfid] = gfid;
//        GlobalFace2LocalFace[gfid] = lfid;
//        lfid++;
//    }

    // ================================== Faces on Rank =========================================
    //NlocElem             = loc_r_elem.size()+Nel_extra+itel;
    //LocalElem2GlobalVert = new Array<int>(NlocElem,8);
    //LocalElem2LocalVert  = new Array<int>(NlocElem,8);
    //std::vector<double> U0vert;
    //U0Elem               = new Array<double>(NlocElem,1);
    //U0Vert               = new Array<double>(LocalVerts.size(),1);
    //ElemPart             = new Array<int>(NlocElem,1);

    
    return pm;
    
}



TetrahedraMesh* ExtractTetrahedralMesh(Array<int>* part_global,
                                       std::map<int,std::vector<int> > tetras,
                                       std::map<int,std::vector<int> > ief_part_map,
                                       std::map<int,std::vector<int> > ifn_part_map,
                                       std::map<int,std::vector<int> > ife_part_map,
                                       std::map<int,std::vector<int> > if_ref_part_map,
                                       std::set<int> ushell,
                                       MPI_Comm comm)
{
	TetrahedraMesh* tmesh = new TetrahedraMesh;
	int world_size;
	MPI_Comm_size(comm, &world_size);
	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(comm, &world_rank);
	    
    int shfn = 0;
    int shf = 0;
    int i;
    int nTetras = tetras.size();
	int gvid,gfid;
	std::set<int> ufaces;
	
	std::map<int,std::vector<int> >::iterator ite;
	int r0,r1,el0,el1,pos,ra;
	int ref,nbcfaces;
	int lcv = 0;
	std::map<int,std::vector<int> > ref2bcface;

	int lf  = 0;
	std::map<int,int> sharedFaces;

//  int* red_ielement_nlocs   = new int[world_size];
	int* ielement_offsets     = new int[world_size];
	int* ielement_nlocs       = new int[world_size];

	int elTel = 0;
	
//	std::map<int,int> lV2gV_tets;
//	std::map<int,int> gV2lV_tets;
	std::set<int> gvid_set;
	int nLocalVerts = 0;
	
//	std::map<int,int> lF2gF_tets;
//	std::map<int,int> gF2lF_tets;
	std::set<int> gfid_set;
	int nLocalFaces = 0;
	
	for(ite=tetras.begin();ite!=tetras.end();ite++)
	{
		//key = GID, value global node nmber;
		int gEl = ite->first;
		
		for(int j=0;j<4;j++)
		{
			gfid = ief_part_map[gEl][j];
			if(gfid_set.find(gfid)==gfid_set.end())
			{
				gfid_set.insert(gfid);
//				gF2lF_tets[gvid] = lfid;
//				lF2gF_tets[lvid] = gfid;
				nLocalFaces++;
			}
			
			for(int k=0;k<3;k++)
			{
				gvid = ifn_part_map[gfid][k];
				
				if(gvid_set.find(gvid)==gvid_set.end())
				{
					gvid_set.insert(gvid);
//					gV2lV_tets[gvid] = lvid;
//					lV2gV_tets[lvid] = gvid;
					nLocalVerts++;
				}
			}
			
			
			if(ushell.find(gfid)!=ushell.end())
			{
				ref = 13;
			}
			else
			{
				ref = if_ref_part_map[gfid][0];
			}

			if(ufaces.find(gfid)==ufaces.end())
			{
				ufaces.insert(gfid);

				el0    = ife_part_map[gfid][0];
				el1    = ife_part_map[gfid][1];

				if(ref==2)
				{
					r0 = part_global->getVal(el0,0);
					r1 = part_global->getVal(el1,0);
					
//					if(r0==world_rank && r1==world_rank)
//					{
//						lf++;
//					}
					if(r0==world_rank && r1!=world_rank)
					{
						sharedFaces[gfid] = r0;
						shf++;
					}
					if(r0!=world_rank && r1==world_rank)
					{
						sharedFaces[gfid] = r1;
						shf++;
					}
				}
				
				if(ref!=2)
				{
					//ref2bcface[ref].push_back(gfid);
					nbcfaces++;
					lf++;
				}
			}
		}
		elTel++;
	}
	
	int nSharedFaces   = sharedFaces.size();
	DistributedParallelState* distTetra 	  = new DistributedParallelState(nTetras,comm);
	DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
	DistributedParallelState* distLocalVerts  = new DistributedParallelState(nLocalVerts,comm);
	DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,comm);

	int Nt_shFaces               = distSharedFaces->getNel();
	int* shFace_offsets          = distSharedFaces->getOffsets();
	int* shFace_nlocs            = distSharedFaces->getNlocs();
	int* shFacesIDs              = new int[nSharedFaces];
	int* shFaces_RankIDs         = new int[nSharedFaces];

	int iter = 0;
	std::set<int> UniqueSharedVertsOnRank_set;
	std::vector<int> UniqueSharedVertsOnRank;
	std::vector<int> UniqueSharedVertsOnRank_RankID;
	
	std::map<int,int>::iterator itsf;
	int lvrtid = 0;
	int tel = shFace_offsets[world_rank];
	for(itsf=sharedFaces.begin();itsf!=sharedFaces.end();itsf++)
	{
		shFacesIDs[iter] = itsf->first;
		shFaces_RankIDs[iter] = itsf->second;
		gfid = itsf->first;

		for(int q=0;q<3;q++)
		{
			gvid   = ifn_part_map[gfid][q];
			
			if(UniqueSharedVertsOnRank_set.find(gvid)==UniqueSharedVertsOnRank_set.end())
			{
				UniqueSharedVertsOnRank_set.insert(gvid);
				UniqueSharedVertsOnRank.push_back(gvid);
				UniqueSharedVertsOnRank_RankID.push_back(world_rank);
				lvrtid++;
			}
		}
		tel++;
		iter++;
	}

	int nSharedVerts = UniqueSharedVertsOnRank.size();

	DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
	
	int Nt_shVerts               = distSharedVerts->getNel();
	int* shVerts_nlocs           = distSharedVerts->getNlocs();
	int* shVerts_offsets         = distSharedVerts->getOffsets();
	
	int* TotalSharedVerts        = new int[Nt_shVerts];
	int* TotalSharedVerts_RankID = new int[Nt_shVerts];
	int* TotalSharedFaces        = new int[Nt_shFaces];
	int* TotalSharedFaces_RankID = new int[Nt_shFaces];
	
	// Communicate vert map to all ranks.
	MPI_Allgatherv(&UniqueSharedVertsOnRank[0],
				   nSharedVerts,
				   MPI_INT,
				   TotalSharedVerts,
				   shVerts_nlocs,
				   shVerts_offsets,
				   MPI_INT, comm);
	
	MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
				   nSharedVerts,
				   MPI_INT,
				   TotalSharedVerts_RankID,
				   shVerts_nlocs,
				   shVerts_offsets,
				   MPI_INT, comm);
	
	// Communicate face map to all ranks.
	MPI_Allgatherv(shFacesIDs,
				   nSharedFaces,
				   MPI_INT,
				   TotalSharedFaces,
				   shFace_nlocs,
				   shFace_offsets,
				   MPI_INT, comm);
	
	MPI_Allgatherv(shFaces_RankIDs,
				   nSharedFaces,
				   MPI_INT,
				   TotalSharedFaces_RankID,
				   shFace_nlocs,
				   shFace_offsets,
				   MPI_INT, comm);
	
	int tmp;
	std::map<int,int> f2r;
	std::set<int> f2r_s;
	int* owned_faces = new int[world_size];
	for(int u=0;u<world_size;u++)
	{
		owned_faces[u] = 0;
	}
	
	int* NewGlobVertCountPerRank = new int[world_size];
	int* NewGlobFaceCountPerRank = new int[world_size];

	for(int u=0;u<world_size;u++)
	{
		NewGlobVertCountPerRank[u] = 0;
		NewGlobFaceCountPerRank[u] = 0;
	}
		
	for(int i=0;i<Nt_shFaces;i++)
	{
		int key = TotalSharedFaces[i];
		int val = TotalSharedFaces_RankID[i];
		
		if(f2r_s.find(key)==f2r_s.end())
		{
			f2r_s.insert(key);
			f2r[key]=val;
			NewGlobFaceCountPerRank[val]=NewGlobFaceCountPerRank[val]+1;
		}
		else
		{
			tmp = f2r[key];
			if(val<tmp)
			{
				f2r[key]=val;
				owned_faces[val]=owned_faces[val]+1;
			}
			if(val>tmp)
			{
				f2r[key]=tmp;
				owned_faces[tmp]=owned_faces[tmp]+1;
			}
		}
	}
	
	std::map<int,int> v2r;
	std::set<int> v2r_s;
	int* owned_verts = new int[world_size];
	for(int u=0;u<world_size;u++)
	{
		owned_verts[u] = 0;
	}
	for(int i=0;i<Nt_shVerts;i++)
	{
		int key = TotalSharedVerts[i];
		int val = TotalSharedVerts_RankID[i];
		
		if(v2r_s.find(key)==v2r_s.end())
		{
			v2r_s.insert(key);
			v2r[key]=val;
			NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
		}
		else
		{
			tmp = v2r[key];
			if(val<tmp)
			{
				v2r[key]=val;
				owned_verts[val]=owned_verts[val]+1;

			}
			if(val>tmp)
			{
				v2r[key]=tmp;
				owned_verts[tmp]=owned_verts[tmp]+1;
			}
		}
	}

	int* NewGlobVertOffset = new int[world_size];
	int* NewGlobFaceOffset = new int[world_size];

	for(int u=1;u<world_size;u++)
	{
		NewGlobVertOffset[u] = NewGlobVertOffset[u-1]+NewGlobVertCountPerRank[u-1];
		NewGlobFaceOffset[u] = NewGlobFaceOffset[u-1]+NewGlobFaceCountPerRank[u-1];
	}
	 
	std::map<int,std::vector<int> >::iterator itv;

	std::map<int,int> sharedVertsGlobal;
	
	for(int u=0;u<world_size;u++)
	{
		NewGlobVertCountPerRank[u] = 0;
		NewGlobFaceCountPerRank[u] = 0;
	}
	
	int iVshared = distLocalVerts->getNel()-v2r.size();

	std::map<int,int >::iterator itvv;
	std::map<int,int> sharedVmap;
	for(itvv=v2r.begin();itvv!=v2r.end();itvv++)
	{
		sharedVmap[itvv->first] = iVshared;
		iVshared++;
	}
	
	std::map<int,int> sharedFmap;
	int iFshared = distLocalFaces->getNel()-f2r.size();

	for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
	{
		sharedFmap[itvv->first] = iFshared;
		iFshared++;
	}
	
	int nNonSharedVerts 		 = nLocalVerts-owned_verts[world_rank];
	int nNonSharedFaces 		 = nLocalFaces-owned_faces[world_rank];

	int* nNonSharedArray         = new int[world_size];
	int* nNonSharedArrayRed      = new int[world_size];
	int* nNonSharedVertsArrayOff = new int[world_size];

	int* nNonSharedFacesArray    = new int[world_size];
	int* nNonSharedFacesArrayRed = new int[world_size];
	int* nNonSharedFacesArrayOff = new int[world_size];

	for(i=0;i<world_size;i++)
	{
		nNonSharedArray[i]      = 0;
		nNonSharedFacesArray[i] = 0;
		if(i==world_rank)
		{
			nNonSharedArray[i] 		= nNonSharedVerts;
			nNonSharedFacesArray[i] = nNonSharedFaces;
		}
	}
	
	MPI_Allreduce(nNonSharedArray,
				  nNonSharedArrayRed,
				  world_size,
				  MPI_INT, MPI_SUM, comm);
	
	MPI_Allreduce(nNonSharedFacesArray,
				  nNonSharedFacesArrayRed,
				  world_size,
				  MPI_INT, MPI_SUM, comm);
	
	int nonSharedOff 		= 0;
	int nonFacesSharedOff 	= 0;
	for(i=0;i<world_size;i++)
	{
		nNonSharedVertsArrayOff[i] = nonSharedOff;
		nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
		
		nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
		nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
	}
	
	//=================================================================================================
	//=================================================================================================
	//=================================================================================================

	int* ini_nEl      = new int[world_size];
	int* red_ini_nEl  = new int[world_size];
	int* ini_offsetEl = new int[world_size];

	for(int i=0;i<world_size;i++)
	{
		ini_nEl[i]      = 0;
		red_ini_nEl[i]  = 0;
		ini_offsetEl[i] = 0;
		
		if(i==world_rank)
		{
			ini_nEl[i] = nTetras;
		}
	}
	MPI_Allreduce(ini_nEl, red_ini_nEl, world_size, MPI_INT, MPI_SUM, comm);

	int offsetEl = 0;
	
	for(int i=0;i<world_size;i++)
	{
		ini_offsetEl[i] = offsetEl;
		offsetEl        = offsetEl + red_ini_nEl[i];
	}
	
	int size = world_size;
	int optimalSize = int(offsetEl/size) + ( world_rank < offsetEl%size );
	double rat = (double)nTetras/(double)optimalSize;
		
	int NtoRecv = 0;
	int NtoSend = 0;
	
	if(nTetras>optimalSize)
	{
		NtoSend = nTetras-optimalSize;
	}
	if(nTetras<optimalSize)
	{
		NtoRecv = optimalSize-nTetras;
	}

   
	int* toS_red = new int[world_size];
	int* toR_red = new int[world_size];
	int* to_Send_copy = new int[world_size];
	int* to_Recv_copy = new int[world_size];
	int* optiSize = new int[world_size];
	int* optiSize_red = new int[world_size];

	int* toS = new int[world_size];
	int* toR = new int[world_size];
	for(int i=0;i<world_size;i++)
	{
		toS[i] = 0;
		toR[i] = 0;
		optiSize[i] = 0;
		optiSize_red[i] = 0;
		if(i==world_rank)
		{
			optiSize[i] = optimalSize;

			toS[i] = NtoSend;
			toR[i] = NtoRecv;
		}
	}
	MPI_Allreduce(optiSize, optiSize_red, world_size, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(toS, toS_red, world_size, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(toR, toR_red, world_size, MPI_INT, MPI_SUM, comm);

	int sent;
	int sendUpdate;
	
	std::map<int,std::vector<int> > recvRa;
	std::map<int,std::vector<int> > recvNe;
	
	std::map<int,std::vector<int> > sendRa;
	std::map<int,std::vector<int> > sendNe;
	
	for(int i=0;i<world_size;i++)
	{
		to_Send_copy[i] =  toS_red[i];
		to_Recv_copy[i] =  toR_red[i];
		if(toR_red[i]!=0)
		{
			for(int j=0;j<world_size;j++)
			{
				if(toS_red[j]>=toR_red[i])
				{
					sendUpdate = toS_red[j]-toR_red[i];
					sent       = toR_red[i];
				}
				
				if(toS_red[j]<toR_red[i])
				{
					sendUpdate = 0;
					sent       = toS_red[j];
				}
				
				toS_red[j] = sendUpdate;
				toR_red[i] = toR_red[i]-sent;
				
				if(sent>0)
				{
					recvRa[i].push_back(j);
					recvNe[i].push_back(sent);
					
					sendRa[j].push_back(i);
					sendNe[j].push_back(sent);
				}
			}
		}
	}
	
	
	//=================================================================================================
	//================================================================================================
	//=================================================================================================
	
	
	int nTet = 0;
	int elloc   = 0;
	int lbvid   = nNonSharedVertsArrayOff[world_rank];
	int lbfid   = nNonSharedFacesArrayOff[world_rank];
	int ltetvid = 0;
	
	std::set<int> gl_set;
	std::map<int,int> gl_map;
	int lvid2   = 0;
	
	std::set<int> glf_set;
	std::map<int,int> glf_map;
	int lbfids  = 0;
	int lfid2   = 0;
	
	int dd      = 0;
	int Nsend   = 0;
	int nv      = 4;
	
	int u = 0;
	int c = 0;
	
	Array<int>* new_ien_or  = new Array<int>(optimalSize,4);
	Array<int>* new_ien     = new Array<int>(optimalSize,4);
	Array<int>* new_ief     = new Array<int>(optimalSize,4);

	if(sendRa.find(world_rank)!=sendRa.end())
	{
		std::vector<int> toRanks    = sendRa[world_rank];
		std::vector<int> NeltoRanks = sendNe[world_rank];
		
		std::vector<std::vector<int> > elNodeIDs;
		std::vector<std::vector<int> > elNodeOriginalIDs;
		std::vector<std::vector<int> > elFaceIDs;
		
		for(int i=0;i<toRanks.size();i++)
		{
			int Nel = NeltoRanks[i];
			std::vector<int> rowNode(Nel*4);
			std::vector<int> rowFace(Nel*4);
			std::vector<int> rowNodeOriginal(Nel*4);

			elNodeIDs.push_back(rowNodeOriginal);
			elNodeOriginalIDs.push_back(rowNode);
			elFaceIDs.push_back(rowFace);
		}
		
		int cc = 0;
		int sRank = toRanks[0];
		
		int offPrank = 0;
		int cntv     = 0;
		
		int t = 0;
		int nuloc = 0;
		int uloc = 0;
		for(ite=tetras.begin();ite!=tetras.end();ite++)
		{	
			int nelPrank  = NeltoRanks[cc];
			int sRank     = toRanks[cc];
			int gEl       = ite->first;
			int lEl       = ini_offsetEl[world_rank]+u;
			int* ien      = new int[4];
			int* ien_o    = new int[4];
			int* ief      = new int[4];
			int* ief_o    = new int[4];
			
			for(int q=0;q<4;q++)
			{
				gvid = ite->second[q];
				
				if(v2r.find(gvid)!=v2r.end())
				{
					lvid2 = sharedVmap[gvid];
					ien[q] = lvid2;
				}
				else
				{
					if(gl_set.find(gvid)==gl_set.end())
					{
						gl_set.insert(gvid);
						gl_map[gvid] = lbvid;
						ien[q] = lbvid;

						lbvid = lbvid + 1;
					}
					else
					{
						int lbbvid = gl_map[gvid];
						ien[q] = lbbvid;
					}
				}
				ien_o[q]=gvid;
			}

			for(int q=0;q<4;q++)
			{
				gfid = ief_part_map[gEl][q];
				
				if(f2r.find(gfid)!=f2r.end())
				{
					lfid2 = sharedFmap[gfid];
					ief[q] = lfid2;
				}
				else
				{
					if(glf_set.find(gfid)==glf_set.end())
					{
						glf_set.insert(gfid);
						glf_map[gfid] = lbfid;
						ief[q] = lbfid;
						lbfid = lbfid + 1;
					}
					else
					{
						int lbbfid = glf_map[gfid];
						ief[q] = lbbfid;
					}
				}
			}
			
			if(u<to_Send_copy[world_rank])
			{
				if(u<(offPrank+nelPrank))
				{
					elNodeIDs[cc][4*t+0]=ien[0];
					elNodeIDs[cc][4*t+1]=ien[1];
					elNodeIDs[cc][4*t+2]=ien[2];
					elNodeIDs[cc][4*t+3]=ien[3];
					
					elNodeOriginalIDs[cc][4*t+0]=ien_o[0];
					elNodeOriginalIDs[cc][4*t+1]=ien_o[1];
					elNodeOriginalIDs[cc][4*t+2]=ien_o[2];
					elNodeOriginalIDs[cc][4*t+3]=ien_o[3];		
										
					elFaceIDs[cc][4*t+0]=ief[0];
					elFaceIDs[cc][4*t+1]=ief[1];
					elFaceIDs[cc][4*t+2]=ief[2];
					elFaceIDs[cc][4*t+3]=ief[3];
					t=t+1;
				}
				else
				{
					t = 0;
					offPrank=offPrank+nelPrank;
					cc=cc+1;
				}
				nuloc++;
			}
			else
			{
				new_ien->setVal(uloc,0,ien[0]);
				new_ien->setVal(uloc,1,ien[1]);
				new_ien->setVal(uloc,2,ien[2]);
				new_ien->setVal(uloc,3,ien[3]);

				new_ien_or->setVal(uloc,0,ien_o[0]);
				new_ien_or->setVal(uloc,1,ien_o[1]);
				new_ien_or->setVal(uloc,2,ien_o[2]);
				new_ien_or->setVal(uloc,3,ien_o[3]);
											
				new_ief->setVal(uloc,0,ief[0]);
				new_ief->setVal(uloc,1,ief[1]);
				new_ief->setVal(uloc,2,ief[2]);
				new_ief->setVal(uloc,3,ief[3]);

				uloc++;
			}
			
			u++;
		}
		
		int acull = 0;
		for(int i=0;i<toRanks.size();i++)
		{
			int dest = toRanks[i];
			int n_El = NeltoRanks[i]*4;
			std::vector<int> Elvec = elNodeIDs[i];
			std::vector<int> ElFvec = elFaceIDs[i];
			std::vector<int> Elovec = elNodeOriginalIDs[i];


			MPI_Send(&n_El  , 1, MPI_INT, dest, dest, comm);
			MPI_Send(&Elvec[0] , n_El, MPI_INT, dest, dest*10000, comm);
			MPI_Send(&ElFvec[0] , n_El, MPI_INT, dest, dest*50000, comm);
			MPI_Send(&Elovec[0] , n_El, MPI_INT, dest, dest*20000, comm);

			acull = acull + n_El;
		}
	}
	
	if(recvRa.find(world_rank)!=recvRa.end())
	{
		std::vector<int > expFromRank = recvRa[world_rank];
		
		std::map<int,std::vector<int> > collected_NIds;
		std::map<int,std::vector<int> > collected_OriginalNIds;
		std::map<int,std::vector<int> > collected_FIds;
		
		for(int i=0;i<expFromRank.size();i++)
		{
			int origin = expFromRank[i];
			int n_Elr;
			MPI_Recv(&n_Elr,   1, MPI_INT, origin, world_rank, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvNElVec(n_Elr);
			MPI_Recv(&recvNElVec[0],   n_Elr, MPI_INT, origin, world_rank*10000, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvFElVec(n_Elr);
			MPI_Recv(&recvFElVec[0],   n_Elr, MPI_INT, origin, world_rank*50000, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvONElVec(n_Elr);
			MPI_Recv(&recvONElVec[0],   n_Elr, MPI_INT, origin, world_rank*20000, comm, MPI_STATUS_IGNORE);
			
			collected_NIds[origin] 			= recvNElVec;
			collected_OriginalNIds[origin] 	= recvONElVec;
			collected_FIds[origin] 			= recvFElVec;
		}
		
		int el = 0;
		
		for(ite=tetras.begin();ite!=tetras.end();ite++)
		{	
			int gEl       = ite->first;
			int lEl       = ini_offsetEl[world_rank]+u;
			int* ien      = new int[nv];
			int* ien_o    = new int[nv];
			int* ief      = new int[nv];
			int* ief_o    = new int[nv];
			
			for(int q=0;q<4;q++)
			{
				gvid = ite->second[q];
				
				if(v2r.find(gvid)!=v2r.end())
				{
					lvid2 = sharedVmap[gvid];
					ien[q] = lvid2;
				}
				else
				{
					if(gl_set.find(gvid)==gl_set.end())
					{
						gl_set.insert(gvid);
						gl_map[gvid] = lbvid;
						ien[q] = lbvid;
						lbvid  = lbvid + 1;
						
					}
					else
					{
						int lbbvid = gl_map[gvid];
						ien[q] = lbbvid;
					}
				}
				ien_o[q] = gvid;
			}
			
			for(int q=0;q<4;q++)
			{
				gfid = ief_part_map[gEl][q];
				
				if(f2r.find(gfid)!=f2r.end())
				{
					lfid2 = sharedFmap[gfid];
					ief[q] = lfid2;
				}
				else
				{
					if(glf_set.find(gfid)==glf_set.end())
					{
						glf_set.insert(gfid);
						glf_map[gfid] = lbfid;
						ief[q] = lbfid;
						lbfid = lbfid + 1;
					}
					else
					{
						int lbbfid = glf_map[gfid];
						ief[q] = lbbfid;
					}
				}
			}
			
			new_ien->setVal(u,0,ien[0]);
			new_ien->setVal(u,1,ien[1]);
			new_ien->setVal(u,2,ien[2]);
			new_ien->setVal(u,3,ien[3]);
			
			new_ien_or->setVal(u,0,ien_o[0]);
			new_ien_or->setVal(u,1,ien_o[1]);
			new_ien_or->setVal(u,2,ien_o[2]);
			new_ien_or->setVal(u,3,ien_o[3]);
			
			new_ief->setVal(u,0,ief[0]);
			new_ief->setVal(u,1,ief[1]);
			new_ief->setVal(u,2,ief[2]);
			new_ief->setVal(u,3,ief[3]);

			u++;
		}
		
		std::map<int,std::vector<int> >::iterator collit;
		int ntot = nTetras;
		for(collit=collected_NIds.begin();collit!=collected_NIds.end();collit++)
		{
			int nel =  collit->second.size()/4;
			ntot = ntot + nel;
			for(int q=0;q<nel;q++)
			{
				new_ien->setVal(u,0,collit->second[q*4+0]);
				new_ien->setVal(u,1,collit->second[q*4+1]);
				new_ien->setVal(u,2,collit->second[q*4+2]);
				new_ien->setVal(u,3,collit->second[q*4+3]);
				
				new_ien_or->setVal(u,0,collected_OriginalNIds[collit->first][4*q+0]);
				new_ien_or->setVal(u,1,collected_OriginalNIds[collit->first][4*q+1]);
				new_ien_or->setVal(u,2,collected_OriginalNIds[collit->first][4*q+2]);
				new_ien_or->setVal(u,3,collected_OriginalNIds[collit->first][4*q+3]);
				
				new_ief->setVal(u,0,collected_FIds[collit->first][4*q+0]);
				new_ief->setVal(u,1,collected_FIds[collit->first][4*q+1]);
				new_ief->setVal(u,2,collected_FIds[collit->first][4*q+2]);
				new_ief->setVal(u,3,collected_FIds[collit->first][4*q+3]);
				
				u++;
			}	
		}	
	}
	
	// return these guys.
	tmesh->ien_part_tetra = new_ien;
	tmesh->ien_part_hybrid = new_ien_or;
	tmesh->ief_part_tetra = new_ief;
    
    return tmesh;
//	new_ien;
//	new_ien_or;
//	new_ief;
	
}


Array<int>* GetNewGlobalPartitioningTetrahedraMesh(TetrahedraMesh* tmesh, MPI_Comm comm)
{
    int i;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    Array<int>* new_ien     = tmesh->ien_part_tetra;
    Array<int>* new_ief     = tmesh->ief_part_tetra;
    Array<int>* new_ien_or  = tmesh->ien_part_hybrid;
    //std::cout << "WORLD_RANK = " << world_rank << " " << lbvid << std::endl;
    
    int optimalSize = new_ien->getNrow();
    int* ielement_nlocs = new int[world_size];
    int* ielement_offsets = new int[world_size];

    int* red_ielement_nlocs = new int[world_size];

    int* elmdist = new int[world_size+1];

    for(i=0;i<world_size;i++)
    {
        ielement_nlocs[i]     = 0;
        red_ielement_nlocs[i] = 0;

        if(i==world_rank)
        {
            ielement_nlocs[i] = optimalSize;
        }
        else
        {
            ielement_nlocs[i] = 0;
        }
    }
    
    MPI_Allreduce(ielement_nlocs,
                  &red_ielement_nlocs[0],
                  world_size, MPI_INT,
                  MPI_SUM, comm);
    
    int o_ie = 0;
    
    for(i=0;i<world_size;i++)
    {
        ielement_offsets[i] = o_ie;
        ielement_nlocs[i]   = red_ielement_nlocs[i];
        elmdist[i]          = o_ie;
        o_ie                = o_ie+red_ielement_nlocs[i];
//        if(world_rank==0)
//        {
//            std::cout << "elmdist[i]  " << elmdist[i] << std::endl;
//
//        }
    }
    
    elmdist[world_size] = o_ie;

    int* eptr     = new int[optimalSize+1];
    int* eind     = new int[optimalSize*4];

    eptr[0]  = 0;
    for(int i=0;i<optimalSize;i++)
    {
        eptr[i+1] = eptr[i]+4;
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j]    = new_ien->data[j];
        }
    }

    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {3};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    int edgecut      = 0;
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];
    idx_t wgtflag_[] = {2};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.02};
    real_t *ubvec    = ubvec_;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    int* part_arr = new int[optimalSize];
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
    real_t itr_[]    = {1.05};
    real_t *itr = itr_;
    
    
    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    int *elmwgt = new int[optimalSize];
    for(int i=0;i<optimalSize;i++)
    {
        elmwgt[i] = 1;
    }

    
    ParMETIS_V3_Mesh2Dual(elmdist,
                          eptr,
                          eind,
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);
    
    ParMETIS_V3_PartKway(elmdist,
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);

    Array<int>* part_global_new  = new Array<int>(o_ie,1);
    Array<int>* part_new         = new Array<int>(optimalSize,1);

    part_new->data = part_arr;
    
    MPI_Allgatherv(&part_new->data[0],
                   red_ielement_nlocs[world_rank], MPI_INT,
                   &part_global_new->data[0],
                   ielement_nlocs,
                   ielement_offsets,
                   MPI_INT,comm);
    
    return part_new;
}

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    
    int ier,opt;
    
    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    const char* fn_metric = "metric.inp";
    
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
    int ReadFromStats = 0;
    if(metric_inputs.size()==6)
    {
        ReadFromStats=metric_inputs[5];
    }
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve = us3d->xcn->getNglob();
    
    int Nel_part = us3d->ien->getNrow();
    
    Array<double>* xcn_ref = ReadDataSetFromFile<double>(fn_grid,"xcn");
    Array<int>* ien_ref = ReadDataSetFromFile<int>(fn_conn,"ien");

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    for(int i=0;i<Nel_part;i++)
    {
        rhoState = us3d->interior->getVal(i,0);
        uState   = us3d->interior->getVal(i,1);
        vState   = us3d->interior->getVal(i,2);
        wState   = us3d->interior->getVal(i,3);
        TState   = us3d->interior->getVal(i,4);
        VtotState = sqrt(uState*uState+vState*vState+wState*wState);
        aState   = sqrt(1.4*287.05*TState);
        MState = VtotState/aState;
        Ui->setVal(i,0,MState);
    }
    
    delete us3d->interior;

    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
      
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);

    
    //std::cout << "understand bitshift " << world_rank << "  " << (world_rank & (~(1 << 3)))<< std::endl;
    //std::vector<int> LocElem     = P->getLocElem();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Domain* pDom = P->getPartitionDomain();
    std::vector<Vert*> VertsOr = P->getLocalVerts();
    
    
    std::map<int,std::vector<int> > tetrasLoc     = pDom->Tetras;
    std::map<int,std::vector<int> > tetras        = pDom->GTetras;
    std::map<int,std::vector<int> > Gtetras       = pDom->GTetras;

    int nTetras       = tetras.size();
    int* Elplease     = new int[world_size];
    int* red_Elplease = new int[world_size];
    std::vector<int> NoElRankID;
    
    for(int u=0;u<world_size;u++)
    {
        Elplease[u]     = 0;
        red_Elplease[u] = 0;
    }
    if(tetras.size()==0)
    {
        Elplease[world_rank] =  1;
    }
    else
    {
        Elplease[world_rank] =  0;
    }
    
    MPI_Allreduce(Elplease, red_Elplease, world_size, MPI_INT, MPI_SUM, comm);
    
    int NRankNoTets = 0;
    int rankalloc   = 0;
    for(int u=0;u<world_size;u++)
    {
        NRankNoTets = NRankNoTets + red_Elplease[u];
        if(red_Elplease[u]==1)
        {
            NoElRankID.push_back(u);
        }
    }
    
    
    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::vector<Vert*> Verts        = P->getLocalVerts();
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    
    std::ofstream myfilet;
    myfilet.open("output_" + std::to_string(world_rank) + ".dat");
    myfilet << "TITLE=\"new_volume.tec\"" << std::endl;
    myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetrasLoc.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid  = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfilet << VertsOr[loc_vid]->x << " "
                << VertsOr[loc_vid]->y << " "
                << VertsOr[loc_vid]->z << std::endl;
    }
    std::map<int,std::vector<int> >::iterator itte;
    for(itte=tetrasLoc.begin();itte!=tetrasLoc.end();itte++)
    {
        myfilet << itte->second[0]+1 << " " << itte->second[1]+1 << " "
                << itte->second[2]+1 << " " << itte->second[3]+1 <<  std::endl;
    }

    myfilet.close();
    
    
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //====================================================================================
  
    //std::vector<int> tris_ref   = pDom->faces_ref;

//    int nTet= 0;
//    std::vector<int> tetcnt = P->getTetCnt();
//    std::map<int,int> mapTets;
//    for(int i=0;i<LocElem.size();i++)
//    {
//        if(tetcnt[i]!=-1)
//        {
//            int gEid = LocElem[i];
//            mapTets[gEid] = tetcnt[i];
//            nTet++;
//        }
//     }
     
    //std::cout << "after partitioning rank = " << world_rank << " #tets = " << tetras.size() << " #prisms " << prisms.size() << std::endl;
    
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* ien_part_map                   = P->getIENpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    Array<int>* part_global                    = P->getGlobalPartition();
    std::set<int> ushell                       = pDom->ushell;
    // map local tets to local tets for parmmg;
    int Nel_glob = part_global->getNrow();
    
    
    TetrahedraMesh* tmesh = ExtractTetrahedralMesh(part_global,tetras,
                                                   ief_part_map->i_map,
                                                   ifn_part_map->i_map,
                                                   ife_part_map->i_map,
                                                   if_ref_part_map->i_map,
                                                   ushell,
                                                   comm);
    
    
    
    
    Array<int>* part_new = GetNewGlobalPartitioningTetrahedraMesh(tmesh,comm);

    
    
    PartMesh* pm = CommunicateTetrahedra(Nel_glob, part_new,
                                         tmesh->ien_part_tetra,
                                         tmesh->ien_part_hybrid,
                                         tmesh->ief_part_tetra,
                                         us3d->xcn, xcn_pstate,
                                         us3d->ifn, ife_pstate, comm);
    
    std::vector<int> lverts;
    std::map<int,int> lpartv2gv_v2;
    std::map<int,int> gv2lpv2;

    std::set<int> gv_set;
    int lcv2 = 0;
    std::map<int,std::vector<int> > ien_part = pm->ien_part;
    std::map<int,std::vector<int> > ien_part_glob = pm->ien_part_glob;
    std::map<int,std::vector<int> >::iterator itertet;
    Array<int>* locelem2locnode= new Array<int>(ien_part.size(),4);
    int elid = 0;
    for(itertet=ien_part.begin();itertet!=ien_part.end();itertet++)
    {
        int glob_id  = itertet->first;
        for(int q=0;q<itertet->second.size();q++)
        {
            int gv = itertet->second[q];
            int lv = pm->GlobalVert2LocalVert[gv];
            
            if(gv_set.find(gv)==gv_set.end())
            {
                gv_set.insert(gv);
                lverts.push_back(lv);
                lpartv2gv_v2[lv]=gv;
                gv2lpv2[gv]=lcv2;
                locelem2locnode->setVal(elid,q,lcv2);
                lcv2=lcv2+1;
            }
            else
            {
                int lcv_u = gv2lpv2[gv];
                locelem2locnode->setVal(elid,q,lcv_u);
            }
        }
        elid++;

    }
    
    
    Array<double>* xcn_glob = ReadDataSetFromFile<double>(fn_grid,"xcn");
    
    std::vector<Vert*> lv = pm->LocalVerts;
    std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << xcn_glob->getNrow() << ", E = " << ien_part.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    std::cout << "loc_part_verts " << lverts.size() << " " << lv.size() << std::endl;
    for(int i=0;i<xcn_glob->getNrow();i++)
    {
        //int loc_vid  = lverts[i];
        //int glob_vid = lpartv2gv_v2[loc_vid];
        //myfile << Verts[loc_vid]->x << " " << Verts[loc_vid]->y << " " << Verts[loc_vid]->z << std::endl;
        myfile << xcn_glob->getVal(i,0) << " " << xcn_glob->getVal(i,1) << " " << xcn_glob->getVal(i,2) << std::endl;
    }
    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
    std::map<int,std::vector<int> >::iterator itt;
    for(itt = ien_part_glob.begin(); itt != ien_part_glob.end(); itt++)
    {
        int glob_id = itertet->first;
        myfile <<   itt->second[0]+1 << "  " <<
        itt->second[1]+1 << "  " <<
        itt->second[2]+1 << "  " <<
        itt->second[3]+1 << "  " << std::endl;
    }


    myfile.close();
    
    
    
    Array<int>* locelem2locnode2= new Array<int>(tmesh->ien_part_hybrid->getNrow(),4);
    std::map<int,int> gv2lpv3;
    std::set<int> gv_setnew;
    elid = 0;
    std::vector<Vert*> locVrts;
    lcv2 = 0;
    for(int i=0;i<tmesh->ien_part_hybrid->getNrow();i++)
	{
		for(int j=0;j<4;j++)
		{
			int gv = tmesh->ien_part_hybrid->getVal(i,j);
			
			if(gv_setnew.find(gv)==gv_setnew.end())
			{
				gv_setnew.insert(gv);
				gv2lpv3[gv]=lcv2;

				Vert* V = new Vert;
				V->x = xcn_glob->getVal(gv,0);
				V->y = xcn_glob->getVal(gv,1);
				V->z = xcn_glob->getVal(gv,2);
				locelem2locnode2->setVal(elid,j,lcv2);
				locVrts.push_back(V);
				
				lcv2=lcv2+1;
			}
			else
			{
				int lcv_u = gv2lpv3[gv];
				locelem2locnode2->setVal(elid,j,lcv_u);
			}
			//locelem2locnode2->setVal(elid,j,gv);

		}
		elid++;

	}
//        
//        
//        
	std::string filename2 = "checkPart_v2_" + std::to_string(world_rank) + ".dat";
	std::ofstream myfile2;
	myfile2.open(filename2);
	myfile2 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
	myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
	myfile2 <<"ZONE N = " << locVrts.size() << ", E = " << tmesh->ien_part_hybrid->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

	for(int i=0;i<locVrts.size();i++)
	{
		myfile2 << locVrts[i]->x << " " << locVrts[i]->y << " " << locVrts[i]->z << std::endl;
	}

	for(int i=0;i<locelem2locnode2->getNrow();i++)
	{
		
		myfile2 <<   locelem2locnode2->getVal(i,0)+1 << "  " <<
					 locelem2locnode2->getVal(i,1)+1 << "  " <<
					 locelem2locnode2->getVal(i,2)+1 << "  " <<
					 locelem2locnode2->getVal(i,3)+1 << "  " << std::endl;
	}


	myfile2.close();
	
	
    /**/
    
    
    
    /**/
    /**/
    
    
//    if(world_rank == 3)
//    {
//        std::map<int,int >::iterator itsf;
//        for(itsf=sharedFacesGlobal.begin();itsf!=sharedFacesGlobal.end();itsf++)
//        {
//            std::cout << itsf->first << " " << itsf->second << std::endl;
//        }
//    }
    
    //std::map<int,std::vector<int> >::iterator itsf;
    
//    for(itsf=recv_SharedFaces.begin();itsf!=recv_SharedFaces.end();itsf++)
//    {
//        std::cout << "recv sizing " << world_rank << " " << itsf->first << " " << itsf->second.size() << std::endl;
//    }
    //std::cout << "recv sizing " << world_rank << " " << recv_SharedFaces.size() << " " << recv_SharedFaces.size() << std::endl;
    
//    ScheduleObj* sObj = DoScheduling(send_SharedFaces,comm);
//
//    std::map<int, std::set<int> >::iterator its;
    
//    if(world_rank == 0)
//    {
//        std::map<int,std::vector<int> >::iterator itm;
//        for(itm=send_SharedFaces.begin();itm!=send_SharedFaces.end();itm++)
//        {
//            std::cout << "Rank  " << world_rank << " sends the following elements: ";
////            for(int q=0;q<itm->second.size();q++)
////            {
////                std::cout <<itm->second[q] << " ";
////            }
//            std::cout << " to rank " << itm->first << std::endl;
//        }
//        //std::cout << "Send Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;
//
//        for(its=sObj->SendFromRank2Rank.begin();its!=sObj->SendFromRank2Rank.end();its++)
//        {
//            std::cout << world_rank << " -> " << its->first << " gets " << its->second.size() << " elements ";
//            std::set<int>::iterator it;
////            for(it=its->second.begin();it!=its->second.end();it++)
////            {
////                std::cout << *it << " ";
////            }
////            std::cout << std::endl;
//        }
////
////
//        //std::cout << "Recv Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;
//
//        for(its=sObj->RecvRankFromRank.begin();its!=sObj->RecvRankFromRank.end();its++)
//        {
//            std::cout << "Rank  " << world_rank << " recieved the following elements: ";
////            for(int q=0;q<itm->second.size();q++)
////            {
////                std::cout <<itm->second[q] << " ";
////            }
//            std::cout << " from rank " << itm->first << std::endl;
//        }
//    }
//
//
//    std::cout << "======================================================="<<std::endl;
//    for(itsf=send_SharedFaces.begin();itsf!=send_SharedFaces.end();itsf++)
//    {
//        std::cout << "send sizing " << world_rank << " " << itsf->first << " " << itsf->second.size() << std::endl;
//    }
    //std::cout << "send sizing " << world_rank << " " << recv_SharedFaces.size() << " " << recv_SharedFaces.size() << std::endl;
    
    
    
//    int* iface_nlocs       = new int[world_size];
//    int* red_iface_nlocs   = new int[world_size];
//    int* iface_offsets     = new int[world_size];
//
//    int* uface_nlocs       = new int[world_size];
//    int* red_uface_nlocs   = new int[world_size];
//    int* uface_offsets     = new int[world_size];
//
//    for(i=0;i<world_size;i++)
//    {
//        iface_nlocs[i]     = 0;
//        red_iface_nlocs[i] = 0;
//
//        uface_nlocs[i]     = 0;
//        red_uface_nlocs[i] = 0;
//
//        if(i==world_rank)
//        {
//            iface_nlocs[i] = uInterFaces.size();
//            uface_nlocs[i] = ufaceOnRank.size();
//        }
//        else
//        {
//            iface_nlocs[i] = 0;
//            uface_nlocs[i] = 0;
//
//        }
//    }
//
//    MPI_Allreduce(iface_nlocs, red_iface_nlocs, world_size, MPI_INT, MPI_SUM, comm);
//    MPI_Allreduce(uface_nlocs, red_uface_nlocs, world_size, MPI_INT, MPI_SUM, comm);
//
//    int o_if = 0;
//    int o_uf = 0;
//
//    for(i=0;i<world_size;i++)
//    {
//        iface_offsets[i] = o_if;
//        o_if = o_if+red_iface_nlocs[i];
//
//        uface_offsets[i] = o_uf;
//        o_uf = o_uf+red_uface_nlocs[i];
//
//    }
//
//    int n_glob_if = o_if;
//
//    std::vector<int> intFace(n_glob_if);
//
//    MPI_Allgatherv(&uSharedFaces[0],
//                   uSharedFaces.size(),
//                   MPI_INT,
//                   &intFace[0],
//                   red_iface_nlocs,
//                   iface_offsets,
//                   MPI_INT, comm);
//
//    std::map<int,int> interFaceMap;
//    int val = 0;
//    for(int u=0;u<intFace.size();u++)
//    {
//        int key = intFace[u];
//        if(interFaceMap.find(key)==interFaceMap.end())
//        {
//            interFaceMap[key] = val;
//            val++;
//        }
//    }
    /* */
    // Communicate the face_color_map;
    //MeshTransfer* TrObj = GetMeshTransfer(tetras,ufaces_vec,uverts_vec,comm,info);
    
//    for(int q=0;q<TrObj->UniqueFacesVec.size();q++)
//    {
//        //std::cout << q << " " << TrObj->UniqueFacesVec[q]<<std::endl;
//    }
//    std::map<int,int> glob2loc_PartInterF = TrObj->iFaces->glob2loc_PartInterEntity;
//    std::set<int> uniqueFaces = TrObj->iFaces->UniqueFaces;
//    tetras     = pDom->GTetras;
//    int uf = 0;
//    ufaces.clear();
//    std::set<int> ufaces2;
//    int duppe3 = 0;
//    for(ite=tetras.begin();ite!=tetras.end();ite++)
//    {
//        int gEl = ite->first;
//
//        for(int j=0;j<4;j++)
//        {
//            gfid = ief_part_map->i_map[gEl][j];
//
//            if(glob2loc_PartInterF.find(gfid)!=glob2loc_PartInterF.end())
//            {
//                duppe3++;
//                ufaces2.insert(gfid);
//            }
//
//        }
//    }
//
//
//    if(world_rank == 0)
//    {
////        std::map<int,int>::iterator its;
////        int dupl = 0;
////        for(its=glob2loc_PartInterF.begin();its!=glob2loc_PartInterF.end();its++)
////        {
////            std::cout << its->first << " " << its->second << std::endl;
////        }
////
////        std::set<int>::iterator itss;
////        int ig = 0;
////        for(itss=ufaces_comp.begin();itss!=ufaces_comp.end();itss++)
////        {
////            std::cout << world_rank << " " << ig << " " << *itss << std::endl;
////            ig++;
////        }
//
//    }
//
//
//
//
//    std::cout << "uf = " << world_rank <<  " " << ufaces.size() << " " << glob2loc_PartInterF.size() << " " << ufaces_VV.size() << " " << ufaces2.size() << " " << duppe3 << " " << ufaces2.size() << " dupl = " << " " << ufaces_comp.size() << " " << uniqueFaces.size() << std::endl;
//
//
    
    
//    std::map<int,int> glob2loc_PartInterV = TrObj->iVerts->glob2loc_PartInterEntity;
//
//    std::cout << "rankert = " << world_rank << " " << duppe << " " << glob2loc_PartInterF.size() << std::endl;
    
    //Array<int>* xcn_tet = new Array<int>(uverts_vec.size(),3);
//    for(ite=tetras.begin();ite!=uverts_vec.end();ite++)
//    {
//        int vid = uverts_vec[i];
//
//        xcn_tet->setVal(i,0,Verts[vid]->x);
//        xcn_tet->setVal(i,1,Verts[vid]->y);
//        xcn_tet->setVal(i,2,Verts[vid]->z);
//    }
    
    
//    Array<int>* ien_tet = new Array<int>(tetras.size(),4);
//    Array<int>* iee_tet = new Array<int>(tetras.size(),6);
//    Array<int>* ief_tet = new Array<int>(tetras.size(),4);
//    ufaces.clear();
//    uverts.clear();
//    int lnfid,lnvid;
//    i = 0;
//    for(ite=tetras.begin();ite!=tetras.end();ite++)
//    {
//        int gEl = ite->first;
//
//        for(int j=0;j<4;j++)
//        {
//            gfid = ief_part_map->i_map[gEl][j];
//            lnfid = glob2loc_PartInterF[gfid];
//            ief_tet->setVal(i,j,lnfid);
//
//            gvid = ien_part_map->i_map[gEl][j];
//            lnvid = glob2loc_PartInterV[gfid];
//            ien_tet->setVal(i,j,lnvid);
//        }
//
//        i++;
//    }
    
    
    
//
    
//+
    
    
    //    Array<int>* ien_tet = GatherTetrahedraOnRoot(tetras,comm,info);

    
    
    /**/
//    int loc_iftri = 0;
//    std::map<int,std::vector<int> >::iterator itm;
//    std::vector<int> interfaces;
//    for(itm  = face_color_map.begin();
//        itm != face_color_map.end();
//        itm++)
//    {
//
//        //std::cout << "R = " << world_rank << " " << itm->second.size() << std::endl;
//
//
//        loc_iftri = loc_iftri + itm->second.size();
//
//        for(int u=0;u<itm->second.size();u++)
//        {
//            interfaces.push_back(itm->second[u]);
//        }
//
//    }
//    //std::cout << "rank " << world_rank << " " << loc_iftri << std::endl;
//
//    std::vector<int> duple = FindDuplicates(interfaces);
    //std::cout << "RR " << world_rank << " " << duple.size() << std::endl;
    
    
    
    /*
    
    std::map<int,std::vector<int> >::iterator itb;
    
    int nbcface     = 0;
    int loc_bftri   = 0;
    int loc_bfqua   = 0;
    int bfaceid     = 0;
    int NvpF        = 0;
    
    // Counting boundary triangles.
    int gv0,gv1,gv2;
    std::vector<std::vector<int> > faces_part;
    std::vector<int> faces_ref;
    for(itb=ref2bcface.begin();itb!=ref2bcface.end();itb++)
    {
        ref     = itb->first;

        for(int fid=0;fid<itb->second.size();fid++)
        {
            bfaceid = itb->second[fid];
            NvpF    = if_Nv_part_map->i_map[bfaceid][0];
            gv0 = gv2lpv[ifn_part_map->i_map[bfaceid][0]];
            gv1 = gv2lpv[ifn_part_map->i_map[bfaceid][1]];
            gv2 = gv2lpv[ifn_part_map->i_map[bfaceid][2]];
            
            std::vector<int> face(3);
            face[0] = gv0;
            face[1] = gv1;
            face[2] = gv2;
            faces_part.push_back(face);
            faces_ref.push_back(ref);

            loc_bftri++;
        }
    }
    
    
    
    int* bftris     = new int[world_size];
    int* red_bftris = new int[world_size];
    int* off_bftris = new int[world_size];
    
    
    int* iftris     = new int[world_size];
    int* red_iftris = new int[world_size];
    int* off_iftris = new int[world_size];

    
    int offset_inter;
    int loc_iftri = 0;
    
    // Counting internal partition triangles.

    std::map<int,std::vector<int> >::iterator itm;
    for(itm  = face_color_map.begin();
        itm != face_color_map.end();
        itm++)
    {
        loc_iftri = loc_iftri + itm->second.size();
    }

    for(int i=0;i<world_size;i++)
    {
        bftris[i]=0;
        iftris[i]=0;

        
        if(i==world_rank)
        {
            bftris[i] = loc_bftri;
            iftris[i] = loc_iftri;
        }
    }
    
    
    MPI_Allreduce(iftris,  red_iftris,   world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(bftris,  red_bftris,   world_size, MPI_INT, MPI_SUM, comm);
    
    int offset_btri = 0;
    int offset_itri = 0;
    
    for(int i=0;i<world_size;i++)
    {
        off_bftris[i] = offset_btri;
        off_iftris[i] = offset_itri;
        
        offset_btri   = offset_btri    +  red_bftris[i];
        offset_itri   = offset_itri   +   red_iftris[i];
    }

    int offset_tri_glob = offset_btri;
    int ncomm           = ordered_rank.size();
    int* color_face     = new int[ncomm];
    std::vector<int* > ifc_tria_glob(ncomm);
    std::vector<int* > ifc_tria_loc(ncomm);

    int* ntifc = new int[ncomm];
    std::set<int>::iterator it;
    int u = 0;
    for(itm  = face_color_map.begin();
        itm != face_color_map.end();
        itm++)
    {
        ra              = itm->first;
        pos             = std::distance(ordered_rank.begin(),ordered_rank.find(ra));
        color_face[pos] = ra;
        ntifc[pos]      = itm->second.size();

        int* ifc_tria_glob_pos = new int[ntifc[pos]];
        int* ifc_tria_loc_pos = new int[ntifc[pos]];

//        std::cout << "rank =  " << rank << " ncomm = "
//         << ncomm << " pos = " << pos << " ra = " << ra
//         << " nF = " << itm->second.size() << std::endl;
        
        for(int i=0;i<itm->second.size();i++)
        {
            std::vector<int> face(3);
            face[0] = gv2lpv[ifn_part_map->i_map[itm->second[i]][0]];
            face[1] = gv2lpv[ifn_part_map->i_map[itm->second[i]][1]];
            face[2] = gv2lpv[ifn_part_map->i_map[itm->second[i]][2]];
            faces_part.push_back(face);
            faces_ref.push_back(2);
            ifc_tria_glob_pos[i] = offset_tri_glob+off_iftris[world_rank]+u;
            ifc_tria_loc_pos[i]  = red_bftris[world_rank]+u;
            u++;
        }
        ifc_tria_glob[pos] = ifc_tria_glob_pos;
        ifc_tria_loc[pos] = ifc_tria_loc_pos;
    }
    
    
    
    
    
    PMMG_pParMesh   parmesh;
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);
    */
    
    
    
    

//==================================================================================
//==================================================================================
//==================================================================================
//    int nVertices       = lv2gpv.size();
//    int nTetrahedra     = tetras.size();
//    int nTriangles      = tris.size();
//    int nEdges          = 0;
//    int nPrisms         = 0;
//    int nQuadrilaterals = bfqua;
//
//    std::map<int,int> gv2lpartv = pDom->gv2lpartv;
//    //std::cout << world_rank << " " << bftri << " " << bfqua << std::endl;
//
//    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
//                              nQuadrilaterals,nEdges) != 1 ) {
//      MPI_Finalize();
//      exit(EXIT_FAILURE);
//    }
//    const int nSolsAtVertices = 1; // 3 solutions per vertex
//    int solType[1] = {MMG5_Tensor};
//
//    if ( PMMG_Set_solsAtVerticesSize(parmesh,nSolsAtVertices,nVertices,solType) != 1 ) {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    }
//
//
//    double* tensor = new double[6];
//    for ( k=0; k<nVertices; ++k )
//    {
//          int gvid = lv2gpv[k];
//          int lvid_p = gv2lpartv[gvid];
//
//          double vx = verts[lvid_p]->x;
//          double vy = verts[lvid_p]->y;
//          double vz = verts[lvid_p]->z;
//
//          if ( PMMG_Set_vertex(parmesh,vx,vy,vz, 1.0, k+1) != 1 )
//          {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//          }
//
//
//
//          tensor[0] = hess_vmap[gvid]->getVal(0,0);
//          tensor[1] = hess_vmap[gvid]->getVal(0,1);
//          tensor[2] = hess_vmap[gvid]->getVal(0,2);
//          tensor[3] = hess_vmap[gvid]->getVal(1,1);
//          tensor[4] = hess_vmap[gvid]->getVal(1,2);
//          tensor[5] = hess_vmap[gvid]->getVal(2,2);
//
//
//          if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,1,tensor,k+1) != 1 )
//          {
//             MPI_Finalize();
//             exit(EXIT_FAILURE);
//          }
//
//
//
//    }
//
//
//
//
//    std::map<int,std::vector<int> >::iterator ittet;
//    k = 0;
//    for ( ittet=tetras.begin(); ittet!=tetras.end(); ittet++ )
//    {
//    	if ( PMMG_Set_tetrahedron(parmesh,ittet->second[0]+1,ittet->second[1]+1,
//        		ittet->second[2]+1,ittet->second[3]+1,1.0,k+1) != 1 ) {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//          }
//        k++;
//    }
//
//
//
//
//
//
//
//    std::map<int,std::vector<int> >::iterator itm;
//    k = 0;
//    for(int q=0;q<tris.size();q++)
//    {
//        int vt0 = tris[q][0];
//        int vt1 = tris[q][1];
//        int vt2 = tris[q][2];
//        int ref = tris_ref[q];
//
//        if ( PMMG_Set_triangle(parmesh, vt0+1,vt1+1,vt2+1, ref,k+1) != 1 )
//        {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//        }
//        k++;
//    }
//
//    int API_mode = PMMG_APIDISTRIB_faces;
//
//    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, 0 ) )
//    {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    };
//
//    if( !world_rank )
//    printf("\n--- API mode: Setting face communicators\n");
//
//    /* Set the number of interfaces */
//    //ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
//    int n_face_comm;
//    //ier = PMMG_Set_numberOfFaceCommunicators(parmesh, n_face_comm);
//    int* color_face     = pDom->color_face;
//    std::vector<int* > ifc_tria_loc  = pDom->ifc_tria_loc;
//    std::vector<int* > ifc_tria_glob = pDom->ifc_tria_glob;
//    int* ntifc          = pDom->ntifc;
//    int ncomm           = pDom->ncomm;
//    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
//
//    //std::cout << n_face_comm << std::endl;
//    /* Loop on each interface (proc pair) seen by the current rank) */
//    for(int icomm=0; icomm<ncomm; icomm++ )
//    {
//
//        std::cout << "rank = " << world_rank << " " << color_face[icomm] << " " << ntifc[icomm] << " " << nTriangles << " bftrie = " << bftri << std::endl;
//    	/* Set nb. of entities on interface and rank of the outward proc */
//    	ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
//										   color_face[icomm],
//										   ntifc[icomm]);
//
//
////        // Set local and global index for each entity on the interface
//    	ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
//                                                 ifc_tria_loc[icomm],
//                                                 ifc_tria_glob[icomm], 1 );
//    }
//
//    int niter = 1;
//
//    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) )
//    {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    };
//
//
//    std::cout << "starting the parallel remeshing..." << std::endl;
//    int ierlib = PMMG_parmmglib_distributed( parmesh );
//    std::cout << "finishing the parallel remeshing..." << std::endl;
//
//
    
    MPI_Finalize();
    
}

