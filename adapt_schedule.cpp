#include "adapt_schedule.h"




ScheduleObj* DoScheduling(std::map<int,std::vector<int> > Rank2RequestEntity, MPI_Comm comm)
{
    int i,t;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nRank_RequestEntity               = Rank2RequestEntity.size(); // The number of ranks from which current rank requests vertices/faces/elements.
    int* reduced_nRank_RequestEntity = new int[size]; // Defined memory for a reduced array so that all ranks are going to be aware of which information is required from each other.
    int* arr_nRank_RequestEntity     = new int[size]; // Defining memory to store local requesting information.

    for(i=0;i<size;i++)
    {
        reduced_nRank_RequestEntity[i] = 0;
        
        if(i==rank)
        {
            arr_nRank_RequestEntity[i] = nRank_RequestEntity+1;
        }
        else
        {
            arr_nRank_RequestEntity[i] = 0;
        }
    }

    // This array hold the information so that each rank knows from how many ranks other ranks request entities.
    MPI_Allreduce(arr_nRank_RequestEntity, reduced_nRank_RequestEntity, size, MPI_INT, MPI_SUM, comm);


    int* reduced_nRank_RequestEntity_offset = new int[size];// Define an offset array for the schedule in order to be able to gather the local schedules for each rank to a global schedule array.

    int offset = 0;
    for(i=0;i<size;i++)
    {
        reduced_nRank_RequestEntity_offset[i] = offset;
        offset = offset+reduced_nRank_RequestEntity[i];
    }
    int nTot_RequestEntity = 0;
    int nRank_RequestEntity_p_one = nRank_RequestEntity+1; // This size is added by one since we add the current rank number to the array.

    MPI_Allreduce(&nRank_RequestEntity_p_one, &nTot_RequestEntity, 1, MPI_INT, MPI_SUM, comm);// Determine the total length of the "schedule" array.
    int* sendFromRank2Rank_Entity_Global = new int[nTot_RequestEntity]; // This array is laid out as follows:
    // first the fromRank is noted which is followed by the IDs for the several ranks FromRank is requesting entities.
    int* sendNentityFromRank2Rank_Global = new int[nTot_RequestEntity]; // This array is layout as follows:
    // first the fromRank is noted which is followed by the number of entities is listed for the several ranks FromRank is sending to.

    for(i=0;i<nTot_RequestEntity;i++)
    {
        sendFromRank2Rank_Entity_Global[i]   = 0;
        sendNentityFromRank2Rank_Global[i]   = 0;
    }

    int* ReqRank_fromRank_Entity   = new int[nRank_RequestEntity+1];
    int* ReqNentity_from_rank      = new int[nRank_RequestEntity+1];
    ReqRank_fromRank_Entity[0]     = rank;
    ReqNentity_from_rank[0]        = -1;

    t = 1;
    std::map<int,std::vector<int> >::iterator it2;
    for(it2=Rank2RequestEntity.begin();it2!=Rank2RequestEntity.end();it2++)
    {
        ReqRank_fromRank_Entity[t] = it2->first;
        ReqNentity_from_rank[t]    = it2->second.size();
        t++;
    }

    /* The example could be:
       rank 0 sends element 1 and 2 to rank 1 and element 10 to rank 2
       rank 1 sends element 2 to rank 0 and element 5 to rank 3
       rank 2 send element 4 to rank 0 and element 3 to rank 3
       rank 3 send element 9 to rank 1 and element 5, 7 and 8 to rank 2

       Then sendFromRank2Rank_Entity_Global and sendNentityFromRank2Rank_Global will look like:

       sendFromRank2Rank_Entity_Global[0]  = 0  sendNentityFromRank2Rank_Global[0]  = -1
       sendFromRank2Rank_Entity_Global[1]  = 1  sendNentityFromRank2Rank_Global[0]  =  2
       sendFromRank2Rank_Entity_Global[2]  = 2  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[3]  = 1  sendNentityFromRank2Rank_Global[0]  = -1
       sendFromRank2Rank_Entity_Global[4]  = 0  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[5]  = 3  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[6]  = 2  sendNentityFromRank2Rank_Global[0]  = -1
       sendFromRank2Rank_Entity_Global[7]  = 0  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[8]  = 3  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[9]  = 3  sendNentityFromRank2Rank_Global[0]  = -1
       sendFromRank2Rank_Entity_Global[10] = 1  sendNentityFromRank2Rank_Global[0]  =  1
       sendFromRank2Rank_Entity_Global[11] = 2  sendNentityFromRank2Rank_Global[0]  =  3
    */

    MPI_Allgatherv(&ReqRank_fromRank_Entity[0],
                   nRank_RequestEntity_p_one, MPI_INT,
                   &sendFromRank2Rank_Entity_Global[0],
                   reduced_nRank_RequestEntity, reduced_nRank_RequestEntity_offset,
                   MPI_INT,comm);

    MPI_Allgatherv(&ReqNentity_from_rank[0],
                   nRank_RequestEntity_p_one, MPI_INT,
                   &sendNentityFromRank2Rank_Global[0],
                   reduced_nRank_RequestEntity, reduced_nRank_RequestEntity_offset,
                   MPI_INT,comm);


    ScheduleObj* scheduleObj   = new ScheduleObj;

    //====================================================================
    for(i=0;i<size;i++)
    {
        int of = reduced_nRank_RequestEntity_offset[i];
        int nl = reduced_nRank_RequestEntity[i];
        
        for(int j=of+1;j<of+nl;j++)
        {
            scheduleObj->SendFromRank2Rank[sendFromRank2Rank_Entity_Global[of]].insert(sendFromRank2Rank_Entity_Global[j]);
            
            scheduleObj->RecvRankFromRank[sendFromRank2Rank_Entity_Global[j]].insert(sendFromRank2Rank_Entity_Global[of]);
        }
    }
    //====================================================================
    delete[] reduced_nRank_RequestEntity_offset;
    delete[] sendNentityFromRank2Rank_Global;
    delete[] sendFromRank2Rank_Entity_Global;
    delete[] ReqRank_fromRank_Entity;
    delete[] ReqNentity_from_rank;
    delete[] reduced_nRank_RequestEntity;
    delete[] arr_nRank_RequestEntity;

    return scheduleObj;
    
}



//void GetAdjacentElementForRank(int *xadj, int* adjcny)
//{
//    std::map<int,std::vector<int> > req_elem;
//
//    for(int i=0;i<part->getNrow();i++)
//    {
//        int start = xadj[i];
//        int end   = xadj[i+1];
//        for(int j=start;j<end;j++)
//        {
//            int adjEl_id = adjcny[j];
//
//            if(elem_set.find(adjEl_id)==elem_set.end())
//            {
//                elem_set.insert(adjEl_id);
//                p_id = part_global->getVal(adjEl_id,0);
//                req_elem[p_id].push_back(adjEl_id);
//            }
//        }
//    }
//
//    std::map<int,std::vector<int> >::iterator itv;
//    for(itv=req_elem.begin();itv!=req_elem.end();itv++)
//    {
//        std::cout <<"rank " << rank << " requests " << " from rank " << itv->first << " the elements -> ";
//        for(int j=0;j<itv->second.size();j++)
//        {
//            std::cout << itv->second[j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//    ScheduleObj* sobj_el = DoScheduling(req_elem,comm);
//
//    std::map<int,std::vector<int> >  reqstd_adj_ids_per_rank;
//    int n_reqstd_adj_ids;
//    for(q=0;q<size;q++)
//    {
//        if(rank==q)
//        {
//            int i=0;
//            for (it = req_elem.begin(); it != req_elem.end(); it++)
//            {
//                int n_req_adj_el           = it->second.size();
//                int dest                   = it->first;
//
//                int destination = dest;
//                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
//                MPI_Send(&n_req_adj_el, 1, MPI_INT, dest, 9876000+10*dest, comm);
//                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
//                MPI_Send(&it->second[0], n_req_adj_el, MPI_INT, dest, 9876000*2+dest*2, comm);
//
//                i++;
//            }
//        }
//        else if (sobj_el->SendFromRank2Rank[q].find( rank ) != sobj_el->SendFromRank2Rank[q].end())
//        {
//            MPI_Recv(&n_reqstd_adj_ids, 1, MPI_INT, q, 9876000+10*rank, comm, MPI_STATUS_IGNORE);
//            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
//
//            std::vector<int> recv_reqstd_adj_ids(n_reqstd_adj_ids);
//            MPI_Recv(&recv_reqstd_adj_ids[0], n_reqstd_adj_ids, MPI_INT, q, 9876000*2+rank*2, comm, MPI_STATUS_IGNORE);
//            reqstd_adj_ids_per_rank[q] = recv_reqstd_adj_ids;
//        }
//    }
//
//    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)
//    {
//        std::cout <<"rank " << rank << " received from " << itv->first << " the request for elements -> ";
//        for(int j=0;j<itv->second.size();j++)
//        {
//            std::cout << itv->second[j] << " ";
//        }
//        std::cout << std::endl;
//    }
//
//    int offset_adj_xcn = 0;
//    int nloc_adj_xcn = 0;
//    std::map<int,int  > recv_adj_back_Nverts;
//    std::map<int,int* > recv_adj_back_verts;
//    std::map<int,int* > recv_adj_back_verts_ids;
//    int n_adj_recv_back;
//
//    // This sends the right vertices of the requested elements to correct processor.
//    for(q=0;q<size;q++)
//    {
//        if(rank == q)
//        {
//            for (it = reqstd_adj_ids_per_rank.begin(); it != reqstd_adj_ids_per_rank.end(); it++)
//            {
//                int nv_adj_send       = it->second.size();
//                double* vert_adj_send = new double[nv_adj_send*3];
//                offset_adj_xcn        = xcn->getOffset(rank);
//                nloc_adj_xcn          = xcn->getNloc(rank);
//                for(int u=0;u<it->second.size();u++)
//                {
//                    vert_adj_send[u*3+0]=xcn->getVal(it->second[u]-offset_adj_xcn,0);
//                    vert_adj_send[u*3+1]=xcn->getVal(it->second[u]-offset_adj_xcn,1);
//                    vert_adj_send[u*3+2]=xcn->getVal(it->second[u]-offset_adj_xcn,2);
//                }
//
//                int dest = it->first;
//                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 98760000+1000*dest, comm);
//                MPI_Send(&it->second[0], it->second.size(), MPI_INT, dest, 9999*9876+dest*8888,comm);
//                delete[] vert_adj_send;
//            }
//        }
//        if(sobj_el->RecvRankFromRank[q].find( rank ) != sobj_el->RecvRankFromRank[q].end())
//        {
//            MPI_Recv(&n_adj_recv_back, 1, MPI_INT, q, 98760000+1000*rank, comm, MPI_STATUS_IGNORE);
//
//            int* recv_adj_back_arr_ids = new int[n_adj_recv_back];
//            MPI_Recv(&recv_adj_back_arr_ids[0], n_adj_recv_back, MPI_INT, q, 9999*9876+rank*8888, comm, MPI_STATUS_IGNORE);
//
//            recv_adj_back_Nverts[q]     = n_adj_recv_back;
//            recv_adj_back_verts_ids[q]  = recv_adj_back_arr_ids;
//        }
//    }
//}
